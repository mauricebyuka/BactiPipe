"""Conditional Illumina filtering and full-precision selected-read provenance."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import shutil
import tempfile

SCHEMA_VERSION = 1


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)


def file_signature(path):
    path = Path(path).resolve()
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def combined_metrics(quality):
    bases = int(quality.combined_total_bases)
    reads = int(quality.combined_num_reads)
    quality_sum = int(quality.combined_quality_sum)
    if bases <= 0 or reads <= 0:
        raise ValueError('No usable reads remain for Illumina analysis.')
    return dict(total_bases=bases, total_reads=reads, quality_sum=quality_sum,
                mean_quality=quality_sum / bases)


def needs_filtering(mean_quality, minimum_quality):
    if not math.isfinite(mean_quality) or not math.isfinite(minimum_quality):
        raise ValueError('Illumina quality thresholds and measurements must be finite.')
    return mean_quality < minimum_quality + 1.0


def select_reads(inputs, qc_root, sample, minimum_quality=28, threads=4, logfile=None):
    from bactipipe.scripts.qualityProc import ProcQuality
    if not math.isfinite(minimum_quality) or not 0 <= minimum_quality <= 93:
        raise ValueError('Illumina minimum quality must be between 0 and 93.')
    if len(inputs) not in (1, 2):
        raise ValueError('Illumina analysis requires one read file or one paired read set.')
    qc_root = Path(qc_root)
    root = qc_root / sample
    root.mkdir(parents=True, exist_ok=True)
    checkpoint = root / 'illumina_read_selection.json'
    contract = dict(schema_version=SCHEMA_VERSION, inputs=[file_signature(p) for p in inputs],
                    minimum_quality=minimum_quality, trigger_quality=minimum_quality+1,
                    implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    quality_implementation_sha256=hashlib.sha256((Path(__file__).parent / 'scripts' / 'qualityProc.py').read_bytes()).hexdigest())
    if checkpoint.is_file():
        old = json.loads(checkpoint.read_text())
        if old.get('contract') == contract:
            try:
                intact = old['selected_files'] == [file_signature(p) for p in old['selected_paths']]
            except OSError:
                intact = False
            if intact and all((qc_root / f'{sample}{suffix}').is_file() for suffix in
                              ('_metrics.txt', '_metrics.png', '_quality_distribution.png')):
                if not old['filter_applied'] or old['fastp_version'] == fastp_version():
                    return old

    (root / 'illumina_filter' / 'invocation.json').unlink(missing_ok=True)
    quality = ProcQuality(list(inputs), None, str(qc_root / sample), sample, logfile=logfile)
    raw = combined_metrics(quality)
    if len(inputs) == 2 and len({s['num_reads'] for s in quality.individual_stats}) != 1:
        raise ValueError('Illumina R1 and R2 contain different numbers of reads.')
    filtered = needs_filtering(raw['mean_quality'], minimum_quality)
    paths = [str(p) for p in inputs]
    command, version = [], None
    if filtered:
        filtering = root / 'illumina_filter'
        filtering.mkdir(exist_ok=True)
        paths = [str(filtering / f'selected_R{i+1}.fastq.gz') for i in range(len(inputs))]
        # fastp accepts integral per-read average-Q gates. Round upward for a
        # fractional configured minimum; retain the exact sample-level gate.
        command = ['fastp', '-i', str(inputs[0]), '-o', paths[0],
                   '--average_qual', str(math.ceil(minimum_quality)),
                   '--qualified_quality_phred', '15', '--unqualified_percent_limit', '40',
                   '--disable_length_filtering', '--thread', str(max(1, min(16, int(threads)))),
                   '--json', str(filtering / 'fastp.json'), '--html', str(filtering / 'fastp.html')]
        if len(inputs) == 2:
            command += ['-I', str(inputs[1]), '-O', paths[1]]
        version = fastp_version()
        atomic_json(filtering / 'invocation.json', dict(command=command, version=version))
        with (filtering / 'fastp.log').open('w') as log:
            log.write('COMMAND ' + json.dumps(command) + '\n')
            log.flush()
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=7200)
        quality = ProcQuality(paths, None, str(qc_root / sample), sample, logfile=logfile)
    selected = combined_metrics(quality)
    if len(inputs) == 2 and len({s['num_reads'] for s in quality.individual_stats}) != 1:
        raise ValueError('Selected Illumina R1 and R2 contain different numbers of reads.')
    # Existing image/TSV plot caches may describe raw reads or a previous input.
    for suffix in ('_metrics.txt', '_metrics.png', '_quality_distribution.png'):
        (qc_root / f'{sample}{suffix}').unlink(missing_ok=True)
    quality.plot_quality_distribution()
    quality.plot_quality_metrics()
    # Canonical metrics are combined across both mates and preserve precision.
    (qc_root / f'{sample}_metrics.txt').write_text(
        f"Avg Quality: {selected['mean_quality']:.17g}\nTotal reads: {selected['total_reads']:,}\n"
        f"Total bases: {selected['total_bases']:,}\n"
        f"Read source: {'filtered' if filtered else 'raw'}\n")
    result = dict(contract=contract, raw=raw, selected=selected, filter_applied=filtered,
                  selected_paths=paths, selected_files=[file_signature(p) for p in paths],
                  fastp_command=command, fastp_version=version,
                  selection_id=hashlib.sha256(json.dumps(contract, sort_keys=True).encode()).hexdigest(),
                  quality_verdict='Pass' if selected['mean_quality'] >= minimum_quality else 'Fail')
    atomic_json(checkpoint, result)
    return result


def fastp_version():
    result = subprocess.run(['fastp', '--version'], capture_output=True, text=True, check=True, timeout=30)
    return (result.stdout or result.stderr).strip()


def confident_depth_organism(fields):
    """Closest candidates must never impose an organism-specific depth gate."""
    if str(fields.match_percent).strip().upper() in {'', 'N/A', 'NA', 'NONE'}:
        return None
    from bactipipe.genome_sizes import GenomeSizeCatalog
    name = fields.coverage_organism or fields.identified_organism
    resolution = GenomeSizeCatalog.from_package().resolve(name)
    return resolution.canonical_name if resolution.is_known else None


def assemble_selected(sample, selection, assembly_dir, assembler, version, threads, logfile, *, memory_gb=None):
    """Never reuse an assembly made from a different read-selection contract."""
    from bactipipe.scripts.process_data import assemble
    root = Path(assembly_dir)
    final = root / 'genomes' / f'{sample}.fasta'
    marker = root / 'genomes' / f'{sample}.input.json'
    key = dict(selection_id=selection['selection_id'], assembler=assembler, version=version,
               selected_metrics=selection['selected'], threads=int(threads), memory_gb=memory_gb)
    if final.is_file() and marker.is_file():
        old = json.loads(marker.read_text())
        if old.get('inputs') == key and old.get('assembly_sha256') == hashlib.sha256(final.read_bytes()).hexdigest():
            return str(final)
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f'.illumina-{sample}-', dir=root) as tmp:
        try:
            assemble(sample=sample, reads=selection['selected_paths'], assembly_dir=tmp,
                     assembler=assembler, sequencer='Illumina', cpus=threads, logfile=logfile,
                     memory_gb=memory_gb)
        finally:
            # Failed jobs must retain their tool diagnostics before temporary cleanup.
            saved_log = Path(tmp) / sample / f'{sample}.log'
            if saved_log.is_file():
                log_root = root / sample
                log_root.mkdir(exist_ok=True)
                shutil.copyfile(saved_log, log_root / f'{sample}.log')
        candidate = Path(tmp) / 'genomes' / f'{sample}.fasta'
        if not candidate.is_file() or not candidate.stat().st_size:
            raise RuntimeError('Genome assembly did not produce a usable FASTA file.')
        final.parent.mkdir(parents=True, exist_ok=True)
        staging = final.with_suffix('.tmp')
        shutil.copyfile(candidate, staging)
        staging.replace(final)
        atomic_json(marker, dict(inputs=key, assembly_sha256=hashlib.sha256(final.read_bytes()).hexdigest()))
    return str(final)


def refresh_depth_row(row, fields, summary_path):
    """Keep depth gates consistent after cohort review and selective reanalysis.

    The new selected-read contract is an explicit Illumina-only opt-in boundary;
    older reports and Nanopore rows are not reinterpreted by this helper.
    """
    from bactipipe.scripts.qc_plan import illumina_minimum_coverage
    from bactipipe.assembly_depth import assembly_depth
    sample = row.get('Sample', row.get('Sample ID', ''))
    summary_path = Path(summary_path)
    root = summary_path.parent.parent if summary_path.parent.name == 'qc_out' else summary_path.parent
    audit = root / 'qc_out' / sample / 'illumina_read_selection.json'
    if not audit.is_file():
        return
    selected = json.loads(audit.read_text())['selected']
    organism = confident_depth_organism(fields)
    threshold = illumina_minimum_coverage(selected['mean_quality'], organism)
    depth, verdict = assembly_depth(selected['total_bases'], root / 'assemblies' / 'genomes' / f'{sample}.fasta',
                                    threshold, audit_path=audit.parent / 'depth_metrics.json')
    if 'min_cov' in row:
        row.update(min_cov=str(threshold), cov_verdict=verdict)
    if 'Required depth' in row:
        row['Required depth'] = f'{threshold}X'
    row['Coverage Depth'] = f'{depth:.2f}X' if depth is not None else 'N/A'
