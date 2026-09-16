"""Competitive reference-breadth confirmation for bacterial read signals."""

from __future__ import annotations

import csv
import gzip
import os
import subprocess
from pathlib import Path
from typing import Sequence

from bactipipe.species_identification import (
    ReferenceCatalog,
    ReferenceRecord,
    TaxonomicSignal,
)


def confirm_signal_breadth(
    *,
    signals: Sequence[TaxonomicSignal],
    catalog: ReferenceCatalog,
    reads: Sequence[str | Path],
    platform: str,
    output_directory: str | Path,
    threads: int = 1,
    maximum_species: int = 5,
    preferred_species_taxids: Sequence[str] = (),
    minimum_mapping_quality: int = 20,
    minimap2_executable: str = "minimap2",
    samtools_executable: str = "samtools",
) -> tuple[TaxonomicSignal, ...]:
    """Map reads competitively to one canonical reference per candidate species.

    Only primary alignments contribute to breadth.  This is intentionally more
    conservative than mapping independently to each species, which can make the
    same conserved reads appear to support several closely related organisms.
    """

    if len(reads) not in {1, 2}:
        raise ValueError("Breadth confirmation requires one read file or an Illumina pair")
    read_paths = [Path(path) for path in reads]
    for path in read_paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Read file is missing or empty: {path}")
    normalised_platform = platform.strip().casefold()
    if normalised_platform not in {"illumina", "nanopore"}:
        raise ValueError("platform must be Illumina or Nanopore")

    preferred = {str(value) for value in preferred_species_taxids}
    ordered_signals = sorted(
        signals,
        key=lambda signal: (
            signal.species_taxid in preferred,
            signal.nomination_score or 0.0,
            signal.read_count,
            signal.read_percent,
        ),
        reverse=True,
    )
    selected: list[tuple[TaxonomicSignal, ReferenceRecord]] = []
    for signal in ordered_signals:
        reference = catalog.canonical_for_species(signal.species_taxid)
        if reference is not None:
            selected.append((signal, reference))
        if len(selected) >= maximum_species:
            break
    if not selected:
        return tuple(signals)

    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    combined_reference = output / "competitive_species_references.fna"
    contig_taxids = _write_combined_reference(
        combined_reference, catalog=catalog, selected=selected
    )
    bam = output / "competitive_species_alignments.bam"
    _map_primary_alignments(
        combined_reference=combined_reference,
        reads=read_paths,
        platform=normalised_platform,
        bam=bam,
        threads=threads,
        output_directory=output,
        minimap2_executable=minimap2_executable,
        samtools_executable=samtools_executable,
        minimum_mapping_quality=minimum_mapping_quality,
    )
    coverage_path = output / "competitive_species_coverage.tsv"
    _run_coverage(
        bam=bam,
        output_path=coverage_path,
        samtools_executable=samtools_executable,
    )
    metrics = _species_mapping_metrics(coverage_path, contig_taxids)
    total_mapped_reads = sum(value[1] for value in metrics.values())
    return tuple(
        TaxonomicSignal(
            species_taxid=signal.species_taxid,
            species_name=signal.species_name,
            read_count=signal.read_count,
            read_percent=signal.read_percent,
            reference_breadth=(
                metrics[signal.species_taxid][0]
                if signal.species_taxid in metrics
                else None
            ),
            # SpeciesFinder does not report Kraken-style minimizers. After
            # mapping-quality-filtered competitive assignment, mapped reads
            # provide the independent support gate. Kraken retains its stricter
            # distinct-minimizer value when one is available.
            distinct_minimizers=(
                signal.distinct_minimizers
                if signal.distinct_minimizers is not None
                else (
                    metrics[signal.species_taxid][1]
                    if signal.species_taxid in metrics
                    else None
                )
            ),
            competitive_mapped_reads=(
                metrics[signal.species_taxid][1]
                if signal.species_taxid in metrics
                else None
            ),
            competitive_mapped_percent=(
                100.0 * metrics[signal.species_taxid][1] / total_mapped_reads
                if signal.species_taxid in metrics and total_mapped_reads
                else (0.0 if signal.species_taxid in metrics else None)
            ),
            nomination_score=signal.nomination_score,
            nomination_percent=signal.nomination_percent,
        )
        for signal in signals
    )


def _write_combined_reference(
    path: Path,
    *,
    catalog: ReferenceCatalog,
    selected: Sequence[tuple[TaxonomicSignal, ReferenceRecord]],
) -> dict[str, str]:
    temporary = path.with_suffix(path.suffix + ".tmp")
    contig_taxids: dict[str, str] = {}
    with temporary.open("w") as destination:
        for _signal, record in selected:
            source = catalog.fasta_for(record)
            sequence_number = 0
            opener = gzip.open if source.suffix == ".gz" else open
            with opener(source, "rt") as handle:
                for line in handle:
                    if line.startswith(">"):
                        sequence_number += 1
                        original = line[1:].strip().split()[0]
                        contig = f"taxid_{record.species_taxid}__{record.accession}__{sequence_number}"
                        contig_taxids[contig] = record.species_taxid
                        destination.write(f">{contig} source={original}\n")
                    else:
                        destination.write(line)
            if sequence_number == 0:
                raise ValueError(f"Reference FASTA contains no records: {source}")
    os.replace(temporary, path)
    return contig_taxids


def _map_primary_alignments(
    *,
    combined_reference: Path,
    reads: Sequence[Path],
    platform: str,
    bam: Path,
    threads: int,
    output_directory: Path,
    minimap2_executable: str,
    samtools_executable: str,
    minimum_mapping_quality: int,
) -> None:
    preset = "sr" if platform == "illumina" else "map-ont"
    mapper_command = [
        minimap2_executable,
        "-ax",
        preset,
        "-t",
        str(max(1, int(threads))),
        str(combined_reference),
        *[str(path) for path in reads],
    ]
    # 0x4 unmapped, 0x100 secondary, 0x800 supplementary.  Excluding these
    # prevents one read from inflating several candidate species.
    view_command = [
        samtools_executable,
        "view",
        "-b",
        "-q",
        str(max(0, int(minimum_mapping_quality))),
        "-F",
        "2308",
        "-",
    ]
    temporary_bam = bam.with_suffix(bam.suffix + ".tmp")
    sort_command = [
        samtools_executable,
        "sort",
        "-@",
        str(max(1, int(threads))),
        "-o",
        str(temporary_bam),
        "-",
    ]
    stderr_path = output_directory / "competitive_mapping.stderr.log"
    with stderr_path.open("w") as stderr:
        mapper = subprocess.Popen(mapper_command, stdout=subprocess.PIPE, stderr=stderr)
        assert mapper.stdout is not None
        view = subprocess.Popen(
            view_command, stdin=mapper.stdout, stdout=subprocess.PIPE, stderr=stderr
        )
        mapper.stdout.close()
        assert view.stdout is not None
        sorter = subprocess.Popen(sort_command, stdin=view.stdout, stderr=stderr)
        view.stdout.close()
        sorter_status = sorter.wait()
        view_status = view.wait()
        mapper_status = mapper.wait()
    if mapper_status or view_status or sorter_status:
        temporary_bam.unlink(missing_ok=True)
        raise RuntimeError(
            "Competitive species mapping failed; see "
            f"{stderr_path} (minimap2={mapper_status}, view={view_status}, sort={sorter_status})"
        )
    if not temporary_bam.is_file():
        raise RuntimeError("Competitive mapping completed without producing a BAM")
    os.replace(temporary_bam, bam)


def _run_coverage(*, bam: Path, output_path: Path, samtools_executable: str) -> None:
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    with temporary.open("w") as output:
        completed = subprocess.run(
            [samtools_executable, "coverage", str(bam)],
            stdout=output,
            stderr=subprocess.PIPE,
            text=True,
        )
    if completed.returncode != 0:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(
            f"samtools coverage failed with exit {completed.returncode}: "
            f"{(completed.stderr or 'unknown error').strip()}"
        )
    os.replace(temporary, output_path)


def _species_mapping_metrics(
    path: Path, contig_taxids: dict[str, str]
) -> dict[str, tuple[float, int]]:
    totals: dict[str, list[int]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"#rname", "startpos", "endpos", "numreads", "covbases"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "samtools coverage output is missing columns: " + ", ".join(sorted(missing))
            )
        for row in reader:
            contig = row["#rname"]
            taxid = contig_taxids.get(contig)
            if taxid is None:
                raise ValueError(f"Coverage output contains unknown reference contig {contig!r}")
            length = int(row["endpos"]) - int(row["startpos"]) + 1
            covered = int(row["covbases"])
            mapped_reads = int(row["numreads"])
            current = totals.setdefault(taxid, [0, 0, 0])
            current[0] += covered
            current[1] += length
            current[2] += mapped_reads
    return {
        taxid: (100.0 * covered / length if length else 0.0, mapped_reads)
        for taxid, (covered, length, mapped_reads) in totals.items()
    }


def _species_breadth(path: Path, contig_taxids: dict[str, str]) -> dict[str, float]:
    """Backward-compatible breadth-only view used by older callers/tests."""

    return {
        taxid: values[0]
        for taxid, values in _species_mapping_metrics(path, contig_taxids).items()
    }
