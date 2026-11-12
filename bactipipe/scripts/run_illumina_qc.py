import os
import sys
import argparse
import csv
import time
import socket
import shutil
import subprocess
import json
import math
from dataclasses import asdict
import hashlib
import tempfile
from pathlib import Path
from datetime import datetime
from bactipipe.scripts import find_organism
from bactipipe.scripts import process_data
from bactipipe.scripts import sample_failures as sample_failure_audit
from bactipipe.scripts.sample_update import run_qc_sample_update
from bactipipe.illumina_filtering import confident_depth_organism
from bactipipe.scripts.utils import time_print, simple_print, logger, pipeheader, excel_reader, download_s3
from bactipipe.progress import emit
from collections import OrderedDict
from bactipipe.__version__ import __version__
from bactipipe.reproducibility import database_identity, database_identity_label, database_log_record
from bactipipe.genome_sizes import GenomeSizeCatalog
from bactipipe.illumina_inputs import stage_s3_reads
from bactipipe.illumina_parallel import (
    SampleJob, cpu_budget, available_memory_gb, plan_resources, prepare_batch, atomic_json,
)
from bactipipe.assembly_depth import assembly_depth
from bactipipe.lambda_control import is_lambda_control, control_summary_row, assess_control, add_control_arguments
from bactipipe.scripts.qc_plan import illumina_minimum_coverage, write_execution_plan
from bactipipe.species_workflow import (
    command_version,
    identification_database_kind,
    identification_database_version,
    identify_assembly,
    legacy_summary_fields,
    python_package_version,
)
from bactipipe.batch_species import (
    reconcile_batch_identifications,
    update_summary_identifications,
)

# Argument Parsing
class CustomHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs['width'] = 100  # Set the width of the help message
        super().__init__(*args, **kwargs)

    def _get_help_string(self, action):
        if action.default is not None and action.default != argparse.SUPPRESS:
            return super()._get_help_string(action)
        return action.help

parser = argparse.ArgumentParser(
    description="Perform conditional read filtering and quality control on Illumina reads, assemble genomes, and assess assembly-based yield depth.",
    formatter_class=CustomHelpFormatter,
    add_help=False  # Disable default help to add it manually
)

required_args = parser.add_argument_group('Required arguments')
optional_args = parser.add_argument_group('Optional arguments')


required_args.add_argument("-n", "--name",
                    help='Name of the person running the pipeline. Needs quotations to handle spaces.\
                        Example: -n "Jane Doe"', required=True)

required_args.add_argument("-r", "--run_name",
                    help="Run name (Run ID)", required=True)

required_args.add_argument("-d", "--fastq_dir",
                    help="Folder containing subfloders of raw fastq files from Illumina instrument.", required=True)

required_args.add_argument("-s", "--storage_type", choices=["local", "s3", "network_mount"], default="local",
                    help="Type of storage where the raw data is located. Default: local")

required_args.add_argument("-l", "--sample_sheet",
                    help="A two-column file with sample (column 1) and organism species (column2).", required=True)

optional_args.add_argument("-o", "--outdir",
                    help="Path to output directory. Default: current directory")
optional_args.add_argument(
    "--input-run-name",
    help="Physical source run name for S3 reads; defaults to --run_name. Outputs keep --run_name.",
)

optional_args.add_argument ("-a", "--assembler", help="Genome assembler to use.",  choices=["Spades", "Unicycler", "Skesa"], default="Skesa")

optional_args.add_argument("-t", "--threads",
                    type=int, help="Total run CPU budget, capped by affinity/container quota. Default: all available")
optional_args.add_argument('--sample-jobs', type=int, default=3,
                    help='Maximum concurrent Illumina QC/assembly samples (SKESA/SPAdes).')
optional_args.add_argument('--sample-threads', type=int, default=20,
                    help='Maximum assembler threads per sample within the total CPU budget.')
optional_args.add_argument('--memory-gb', type=float,
                    help='Maximum run preparation RAM budget in GiB; automatically bounded by available resources.')
optional_args.add_argument("-q", "--minqual", help="Minimum selected-read mean Phred Q. Filter both raw mates when their combined mean is below this value + 1.", default=28, type=float)
optional_args.add_argument("--min_completeness", help="Minimum CheckM completeness percentage.", default=90, type=float)
optional_args.add_argument("--max_contamination", help="Maximum CheckM contamination percentage.", default=10, type=float)
optional_args.add_argument(
    "--species-database-root",
    help=(
        "Reviewed BactiPipe species-identification database version. The new "
        "ANI workflow remains disabled when omitted."
    ),
)
optional_args.add_argument(
    "--species-gtdb-metadata",
    help="GTDB metadata TSV used to label an unresolved sample's bac120 tree.",
)
optional_args.add_argument(
    "--species-gtdbtk-conda-env",
    default="gtdbtk",
    help="Conda environment containing GTDB-Tk for unresolved samples.",
)
optional_args.add_argument(
    "--allow-candidate-species-database",
    action="store_true",
    help="Permit an inactive species database for controlled benchmarking.",
)
optional_args.add_argument(
    "--only-sample",
    action="append",
    default=[],
    metavar="SAMPLE_ID",
    help="Process only this sample-sheet sample. Repeat for multiple samples.",
)
optional_args.add_argument("--defer-qc-pdf", action="store_true",
                           help="Write QC tables but defer PDF publication to a follow-up workflow.")
optional_args.add_argument(
    "--update-existing",
    action="store_true",
    help=(
        "Atomically replace the selected samples in an existing run output. "
        "Requires --only-sample."
    ),
)

optional_args.add_argument("--s3_bucket", "--s3-bucket", dest="s3_bucket",
                    default=os.getenv("ANALYSIS_S3_BUCKET"),
                    help="S3 bucket for s3 storage_type. Defaults to ANALYSIS_S3_BUCKET.")

optional_args.add_argument("--s3_prefix", "--s3-prefix", dest="s3_prefix",
                    default=os.getenv("ANALYSIS_S3_PREFIX", ""),
                    help="S3 prefix for s3 storage_type. Defaults to ANALYSIS_S3_PREFIX.")

optional_args.add_argument("-h", "--help",
                           action="help",
                           help="Show this help message and exit.")

add_control_arguments(optional_args)
args = parser.parse_args()
if min(args.sample_jobs, args.sample_threads, args.threads if args.threads is not None else 1) < 1:
    parser.error('Thread and sample concurrency limits must be positive')
if args.memory_gb is not None and (not math.isfinite(args.memory_gb) or args.memory_gb < 4):
    parser.error('--memory-gb must be finite and at least 4 GiB')

species_identification_enabled = bool(args.species_database_root)
if args.species_gtdb_metadata and not species_identification_enabled:
    parser.error("--species-gtdb-metadata requires --species-database-root")

# Version of BactiPipe
bactipipe_version = f"v{__version__}"

sample_list = os.path.abspath(args.sample_sheet)
run_name = args.run_name
raw_reads = os.path.abspath(args.fastq_dir)

# raw_reads = os.path.join(raw_reads, run_name)
tech_name = args.name

outdir_parent = os.path.abspath(args.outdir or os.getcwd())
outDir = os.path.join(outdir_parent, run_name)


def identify_sample_with_ani(sample, expected_organism, genome, reads, threads):
    """Run the opt-in ANI workflow with candidate-only Kraken evidence."""
    result = identify_assembly(
        sample_id=sample,
        assembly=genome,
        database_root=args.species_database_root,
        output_directory=os.path.join(outDir, "species_identification", sample),
        expected_organism=expected_organism,
        reads=reads,
        platform="Illumina",
        kraken_database=(
            None
            if identification_database_kind(args.species_database_root)
            == "speciesfinder_kma"
            else os.path.join(args.species_database_root, "kraken2")
        ),
        threads=max(1, int(threads)),
        allow_candidate_database=args.allow_candidate_species_database,
        require_mixed_signal_assessment=True,
    )
    return legacy_summary_fields(result)

if args.update_existing:
    try:
        run_qc_sample_update(
            script_path=__file__,
            argv=sys.argv[1:],
            outdir_parent=outdir_parent,
            run_name=run_name,
            sample_names=args.only_sample,
            species_database_root=args.species_database_root,
            species_platform="Illumina",
            species_gtdb_metadata=args.species_gtdb_metadata,
            species_gtdbtk_conda_env=args.species_gtdbtk_conda_env,
            species_placement_threads=cpu_budget(args.threads),
        )
    except (OSError, RuntimeError, ValueError) as exc:
        sys.stderr.write(f"BactiPipe sample update failed: {exc}\n")
        sys.exit(1)
    sys.exit(0)

def assembly_version():
    assembler = args.assembler.lower()
    if assembler == "spades":
        assembler = "spades.py"
    elif assembler == "unicycler":
        assembler = "unicycler"
    elif assembler == "skesa":
        assembler = "skesa"
    
    vercmd = [assembler, "--version"]
    try:
        process = subprocess.Popen(vercmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        if assembler == "skesa":
            res_line = out.decode('utf-8', errors='replace').strip().split('\n')[-1]
            assembly_ver = res_line.split()[-1]
        else:
            assembly_ver = out.decode('utf-8', errors='replace').strip().split()[-1]
        if not assembly_ver.lower().startswith("v"):
            assembly_ver = "v" + assembly_ver
    except Exception as e:
        assembly_ver = "Unknown"
    return assembly_ver

def get_fastqs(sample, raw_reads):
    fastqs = []
    fastq1 = None
    fastq2 = None
    if not raw_reads or not os.path.isdir(raw_reads):
        return fastqs

    for item in sorted(os.listdir(raw_reads)):
        if item.startswith(sample):
            item_dir = os.path.join(raw_reads, item)
            if not os.path.isdir(item_dir):
                continue
            for file in sorted(os.listdir(item_dir)):
                if "_R1" in file and (file.endswith("fastq") or file.endswith("fastq.gz")):
                    fastq1 = os.path.join(raw_reads, item, file)
                elif "_R2" in file and (file.endswith("fastq") or file.endswith("fastq.gz")):
                    fastq2 = os.path.join(raw_reads, item, file)

    if fastq1 and not fastq2:  # If the data is single end
        fastqs.append(fastq1)
    elif fastq1 and fastq2:
        fastqs.append(fastq1)
        fastqs.append(fastq2)
    return fastqs

# Run the pipeline

if not os.path.exists(outDir):
    os.makedirs(outDir)

cpus = cpu_budget(args.threads)

static = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")

log_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
log = os.path.join(outDir, f"{log_time}~{run_name}_qc.log")

print("")

date = datetime.now().strftime("%Y-%m-%d")

hostname = socket.gethostname()
try:
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        # doesn't even have to be reachable
        s.connect(('10.254.254.254', 1))
        ip_address = s.getsockname()[0]
    except Exception as e:
        ip_address = f"Unable to resolve hostname: {e}"
    finally:
        s.close()
except Exception as e:
    ip_address = f"Unable to resolve hostname: {e}"
    ip_address = "Unable to resolve hostname"

header_info = pipeheader("Illumina", date, tech_name, hostname, ip_address, run_name, sample_list, raw_reads, outDir, cpus)

header = header_info[0]
run_info = header_info[1]
for line in header:
    simple_print(line)
    logger(log, line, mode="simple")
for line in run_info:
    simple_print(line)
    logger(log, line, mode="simple")

# Ensure that the sample list if formatted correctly
if not os.path.exists(sample_list):
    print(f"Sample list file {sample_list} does not exist.")
    logger(log, f"Sample list file {sample_list} does not exist.")
    sys.exit(1)
else:
    if sample_list.endswith('.xlsx') or sample_list.endswith('.xls'):
        try:
            sample_info = excel_reader(sample_list, "illumina")
        except Exception as e:
            print(f"Error reading Excel file {sample_list}: {e}")
            logger(log, f"Error reading Excel file {sample_list}: {e}")
            sys.exit(1)
    else:
        with open(sample_list, 'r') as sampL:
            sample_info = sampL.readlines()
            first_line = sample_info[0].strip().split('\t')
            if len(first_line) != 2:
                print("Sample list file must have two columns: sample, organism.")
                logger(log, "Sample list file must have two columns: sample, organism.")
                sys.exit(1)

if args.only_sample:
    requested_samples = tuple(
        dict.fromkeys(value.strip() for value in args.only_sample if value.strip())
    )
    selected_lines = []
    found_samples = set()
    for line in sample_info:
        if not line.strip() or line.startswith("#"):
            continue
        sample_id = line.strip().split("\t", 1)[0]
        if sample_id == "sample_ID":
            selected_lines.append(line)
        elif sample_id in requested_samples:
            selected_lines.append(line)
            found_samples.add(sample_id)
    missing_samples = [
        sample_id for sample_id in requested_samples if sample_id not in found_samples
    ]
    if missing_samples:
        message = (
            "Requested sample(s) were not found in the sample sheet: "
            + ", ".join(missing_samples)
        )
        print(message)
        logger(log, message)
        sys.exit(1)
    sample_info = selected_lines

# Stage remote reads before checking sample membership; the mount may be stale.
temp_dir = os.path.join(outDir, "temp")
if args.storage_type == "s3":
    if not args.s3_bucket:
        parser.error("S3 storage selected, but no S3 bucket was provided.")
    logger(log, f"Staging S3 reads for source run: {args.input_run_name or run_name}")
    try:
        raw_reads = stage_s3_reads(
            bucket=args.s3_bucket, prefix=args.s3_prefix, run_name=run_name,
            input_run_name=args.input_run_name, target_root=temp_dir, download=download_s3,
        )
    except (OSError, ValueError) as error:
        logger(log, f"Input staging failed: {error}")
        time_print(f"Input staging failed: {error}", "Fail")
        sys.exit(1)
elif args.storage_type == "network_mount":
    copied_reads = os.path.join(temp_dir, os.path.basename(raw_reads))
    shutil.copytree(raw_reads, copied_reads, dirs_exist_ok=True)
    raw_reads = copied_reads

# Resolve the controlled organism and size once per sample. Illumina assembly
# itself remains size-independent; the size is used only for QC depth.
genome_size_catalog = GenomeSizeCatalog.from_package()

# Validate the sample sheet
bad_organisms = []
bad_samples = []
sample_failures = {}
expected_organisms = {}
sample_resolutions = {}
# with open(sample_list, 'r') as sampL:
for line in sample_info:
    if not line.strip() or line.startswith("#"):
        continue  # Skip comment lines
    sample, organism = line.rstrip('\r\n').split('\t')
    if sample == "sample_ID":
        continue
    expected_organisms[sample] = organism
    resolution = genome_size_catalog.resolve(organism)
    sample_resolutions[sample] = resolution
    if resolution.status in {"unsupported", "missing"}:
        reason = (
            f"Organism [{organism.strip()}] is not in the controlled BactiPipe "
            "organism list or was left blank. Continuing in Unknown mode."
        )
        time_print(f"WARNING: Sample {sample}: {reason}", "Fail")
        logger(log, f"WARNING: Sample {sample}: {reason}")
        bad_organisms.append(f"{sample}:{organism}")

    if not get_fastqs(sample, raw_reads) and sample != "sample_ID":
        time_print(f"WARNING: Sample [{sample}] does not have corresponding fastq files", "Fail")
        logger(log, f"WARNING: Sample [{sample}] does not have corresponding fastq files")
        if sample not in bad_samples:
            bad_samples.append(sample)
        sample_failure_audit.record_failure(
            sample_failures,
            sample,
            "input",
            "No corresponding FASTQ files were found.",
        )

if bad_organisms:
    time_print("Unsupported or blank organisms will use Unknown mode.")
    logger(log, "Unsupported or blank organisms will use Unknown mode.")

if bad_samples:
    b_s = ", ".join(bad_samples)
    time_print(f"WARNING: These samples don't have corresponding fastq data. They will be skipped: {b_s}.", "Fail")
    logger(log, f"WARNING: These samples don't have corresponding fastq data. They will be skipped: {b_s}.")

start_time = time.time()
logger(log, "PIPELINE STARTED", "Header")
time_print("PIPELINE STARTED", "Header")

# Count total number of samples
sample_number = sum(
    1
    for line in sample_info
    if line.strip() and not line.startswith("#") and line.strip().split('\t')[0] != "sample_ID"
) - len(bad_samples)
sample_number = max(0, sample_number)

time_print(f"Total number of samples to be processed: {sample_number}\n")
logger(log, f"Total number of samples to be processed: {sample_number}\n")

if sample_number == 0:
    time_print("No valid samples to process. Exiting pipeline.", "Fail")
    logger(log, "No valid samples to process. Exiting pipeline.")
    sys.exit(1)

write_execution_plan(
    os.path.join(outDir, "sample_execution_plan.tsv"),
    [
        (sample, expected_organisms[sample], sample_resolutions[sample])
        for sample in expected_organisms
    ],
    platform="Illumina",
    min_quality=args.minqual,
    min_read_length=None,
    min_coverage="quality_dependent_30_40_50_plus_ecoli_10",
    catalog_checksum=genome_size_catalog.checksum,
    assembly_uses_genome_size=False,
)

qc_out = os.path.join(outDir, "qc_out")
if not os.path.exists(qc_out):
    os.makedirs(qc_out)

#Write the temporary summary file that will be updated after CheckM analysis
temp_qc_summary = os.path.join(qc_out, "temp_qc_summary.tsv")
reported_samples = set()
read_selections = {}
current_genomes = {}
current_assembly_version = assembly_version()
memory_budget = available_memory_gb()
if args.memory_gb is not None:
    memory_budget = min(memory_budget, args.memory_gb)
resources = plan_resources(sample_number, cpus, maximum_jobs=args.sample_jobs,
                           maximum_threads=args.sample_threads, memory_gb=memory_budget,
                           assembler=args.assembler)
preparation_jobs = [
    SampleJob(sample, tuple(get_fastqs(sample, raw_reads)), qc_out,
              os.path.join(outDir, 'assemblies'), args.assembler,
              current_assembly_version, args.minqual, resources.threads_per_sample,
              resources.memory_per_sample_gb)
    for sample in expected_organisms if sample not in bad_samples
]
atomic_json(Path(outDir) / 'illumina_resources.json', asdict(resources))
def report_preparation(message):
    time_print(message)
    logger(log, message)
    sys.stdout.flush()
emit('inputs', 'completed')
emit('read_qc', total=len(preparation_jobs))
emit('assembly', total=len(preparation_jobs))
prepared_samples = prepare_batch(preparation_jobs, resources, report_preparation)
emit('read_qc', 'completed')
emit('assembly', 'completed')
emit('identification', total=len(preparation_jobs))
# Identification starts only after preparation releases its CPU/RAM budget.
from bactipipe.identification_parallel import identify_batch, result_for
identification_results = {}
if species_identification_enabled:
    identification_jobs = {}
    for sample, preparation in prepared_samples.items():
        organism = expected_organisms[sample]
        genome = preparation.get('genome')
        if (sample in bad_samples or is_lambda_control(organism) or preparation.get('failure')
                or not genome or not os.path.isfile(genome) or not os.path.getsize(genome)):
            continue
        selected = preparation.get('selected') or {}
        if selected.get('selected', {}).get('mean_quality', -1) < args.minqual:
            continue
        resolution = sample_resolutions[sample]
        identification_jobs[sample] = (resolution.canonical_name if resolution.is_known else 'unknown',
                                       genome, selected['selected_paths'])
    identification_results = identify_batch(identification_jobs, identify_sample_with_ani,
        cpus=cpus, database_root=args.species_database_root, memory_gb=min(memory_budget, available_memory_gb()),
        report=report_preparation, audit_path=os.path.join(qc_out, 'identification_resources.json'))
with(open(temp_qc_summary , 'w')) as qc_sum:
    writer = csv.writer(qc_sum, dialect='excel-tab')
    writer.writerow(["Sample",  "Mean_quality", "qc_verdict", "Expected organism", "Identified organism", "% Match", "Coverage Depth", "min_cov", "cov_verdict", "tax_confirm"])

    # with open(sample_list, 'r') as sampL:
    for line in sample_info:
        if not line.strip() or line.startswith("#"):
            continue
        sample, organism = line.rstrip('\r\n').split('\t')
        if sample == "sample_ID":
            continue
        if sample in bad_samples:
            continue
        if sample not in identification_results:
            emit('identification', sample=sample)
        fastqs = get_fastqs(sample, raw_reads)
        if not fastqs:
            time_print(f"WARNING: Sample [{sample}] does not have corresponding fastq files", "Fail")
            logger(log, f"WARNING: Sample [{sample}] does not have corresponding fastq files")
            continue
        resolution = sample_resolutions[sample]
        if not resolution.is_known:
            sys_organism = "unknown"
            genome_size = None
            taxID = "N/A"
        else:
            sys_organism = resolution.canonical_name
            genome_size = resolution.genome_size_bp
            taxID = resolution.taxid
        preparation = prepared_samples[sample]
        if preparation['failure']:
            failure = preparation['failure']
            sample_failure_audit.record_failure(sample_failures, sample, failure['stage'], failure['reason'])
            # Preserve filtering invocation provenance/cleanup even if assembly failed.
            if preparation['selected']:
                read_selections[sample] = preparation['selected']
            continue
        try:
            selected = preparation['selected']
            read_selections[sample] = selected
            fastqs = selected['selected_paths']
            avqc = selected['selected']['mean_quality']
            total_bases = selected['selected']['total_bases']
            logger(log, f"Sample {sample}: raw combined Q={selected['raw']['mean_quality']:.17g}; "
                        f"filter trigger <{args.minqual+1:g}; filtered={selected['filter_applied']}; "
                        f"selected Q={avqc:.17g}; selected bases={total_bases}")
        except Exception as e:
            reason = f"Read quality processing failed: {e}"
            time_print(f"Sample {sample}: {reason}", "Fail")
            logger(log, f"Sample {sample}: {reason}")
            sample_failure_audit.record_failure(sample_failures, sample, "qc", reason)
            continue

        taxonomy = [taxID, organism]

        if avqc >= args.minqual:
            qc_verdict = 'Pass'
        else:
            qc_verdict = "Fail"

        # Coverage
        min_cov = illumina_minimum_coverage(avqc, None)
        coverage, cov_verdict = None, 'N/A'

        quality_message = f'Sample {sample}: mean quality score: {qc_verdict} --- {avqc:.2f}'
        simple_print(quality_message, qc_verdict)
        logger(log, quality_message, qc_verdict, mode="simple")
        cov_display = "N/A"
        
        # Preparation already assembled samples that passed read QC. Assess the
        # retained assemblies here; identification can still adjust required depth.
        assembly_dir = os.path.join(outDir, "assemblies")
        if qc_verdict == "Pass":
            genome = preparation['genome']

            if not genome or not os.path.isfile(genome) or os.path.getsize(genome) == 0:
                sample_failure_audit.record_failure(
                    sample_failures,
                    sample,
                    "assembly",
                    "Genome assembly did not produce a usable FASTA file.",
                )
                continue

            coverage, cov_verdict = assembly_depth(
                total_bases, genome, min_cov,
                audit_path=os.path.join(qc_out, sample, "depth_metrics.json"),
            )
            cov_display = f"{coverage:.2f}X" if coverage is not None else "N/A"
            current_genomes[sample] = genome
            depth_message = f'Sample {sample}: assembly-based yield depth: {cov_display}'
            time_print(depth_message)
            logger(log, depth_message, mode="simple")

            if is_lambda_control(organism):
                control_evidence = assess_control(genome, os.path.join(qc_out, sample, "lambda_control"),
                    min_identity=args.lambda_min_identity, min_coverage=args.lambda_min_coverage,
                    report=lambda message: logger(log, message))
                writer.writerow(control_summary_row(sample, avqc, qc_verdict, organism, cov_display, min_cov, cov_verdict, control_evidence))
                if control_evidence["status"] == "Error":
                    sample_failure_audit.record_failure(sample_failures, sample, "control_confirmation", control_evidence["reason"])
                reported_samples.add(sample)
                emit('identification', 'failed' if control_evidence['status'] == 'Error' else 'completed', sample=sample)
                continue

            if species_identification_enabled:
                try:
                    summary_fields = result_for(identification_results, sample)
                except Exception as e:
                    time_print(f"Error at species identification step: {e}", "Fail")
                    logger(log, f"Error at species identification step: {e}")
                    sample_failure_audit.record_failure(
                        sample_failures,
                        sample,
                        "organism_identification",
                        f"Species identification failed: {e}",
                    )
                    continue
                best_org = summary_fields.identified_organism
                best_percent = summary_fields.match_percent
                tax_confirm = summary_fields.taxonomic_confirmation
                identified_org = confident_depth_organism(summary_fields)
            else:
                try:
                    hit, tax_confirm, possibilities = find_organism.find_species_with_kmrf(s_name=sample, lab_species=sys_organism, genome=genome, dataOut=outDir, org_type="bacteria", logfile=log)
                except Exception as e:
                    time_print(f"Error at find organism step: {e}", "Fail")
                    logger(log, f"Error at find organism step: {e}")
                    sample_failure_audit.record_failure(
                        sample_failures,
                        sample,
                        "organism_identification",
                        f"Organism identification failed: {e}",
                    )
                    continue

                best_org = None
                best_percent = None
                best_other_org = None
                identified_org = None
                try:
                    if tax_confirm == "Pass" or tax_confirm == "N/A":
                        best_org = hit.split(" (")[0]
                        best_percent = f"{float(hit.split(' (')[1].strip(')%')):.2f}"
                        identified_org = best_org
                    else:
                        best_other_hit = hit.split(" --- ")[0].split(": ")[1]
                        best_other_org = best_other_hit.split(" (")[0]
                        best_other_percent = best_other_hit.split(" (")[1].strip(")")
                        identified_org = None  # Closest legacy hits are not confirmed organisms.
                except (AttributeError, IndexError, TypeError, ValueError) as e:
                    sample_failure_audit.record_failure(
                        sample_failures,
                        sample,
                        "organism_identification",
                        f"Organism-identification output could not be parsed: {e}",
                    )
                    continue
            # Identification may adjust the existing organism-specific QC requirement,
            # but never the denominator used to report measured post-assembly yield.
            identified_resolution = genome_size_catalog.resolve(identified_org)
            min_cov = illumina_minimum_coverage(avqc, identified_resolution.canonical_name if identified_resolution.is_known else None)
            coverage, cov_verdict = assembly_depth(
                total_bases, genome, min_cov,
                audit_path=os.path.join(qc_out, sample, "depth_metrics.json"),
            )
            cov_display = f"{coverage:.2f}X" if coverage is not None else "N/A"
            depth_message = (f"Sample {sample}: yield depth QC: {cov_verdict} --- "
                             f"{cov_display}; required {min_cov:g}X")
            time_print(depth_message)
            logger(log, depth_message, mode="simple")

            if species_identification_enabled:
                writer.writerow([sample, f'{avqc:.2f}', qc_verdict, organism, best_org, best_percent, cov_display, min_cov, cov_verdict, tax_confirm])
            elif tax_confirm == "Fail":
                writer.writerow([sample, f'{avqc:.2f}', qc_verdict, organism, f'Closest: {best_other_org}', best_other_percent, f'{cov_display}', min_cov, cov_verdict, tax_confirm])
            else:
                writer.writerow([sample, f'{avqc:.2f}', qc_verdict, organism, best_org, best_percent, cov_display, min_cov, cov_verdict, tax_confirm])
            reported_samples.add(sample)
            emit('identification', 'completed', sample=sample)

        else:
            best_org = "N/A"
            best_percent = "N/A"
            tax_confirm = "N/A"

            # A planning estimate is not a final measured depth without an assembly.
            cov_display, cov_verdict = "N/A", "N/A"

            writer.writerow([sample, f'{avqc:.2f}', qc_verdict, organism, best_org, best_percent, cov_display, min_cov, cov_verdict, tax_confirm])
            reported_samples.add(sample)
            emit('identification', 'skipped', sample=sample)

sample_failure_audit.publish_failures(
    outDir,
    temp_qc_summary,
    expected_organisms,
    sample_failures,
    required_depth=0,
)
if not reported_samples:
    time_print("No samples completed analysis. Exiting pipeline.", "Fail")
    logger(log, "No samples completed analysis. Exiting pipeline.")
    sys.exit(1)
# Assess Genome Quality with CheckM
for progress_sample in expected_organisms:
    if progress_sample not in bad_samples and progress_sample not in reported_samples:
        emit('identification', 'failed' if sample_failures.get(progress_sample, ('', ''))[0] == 'organism_identification' else 'skipped', sample=progress_sample)
emit('identification', 'completed')
emit('assembly_qc')
genomes_dir = os.path.join(outDir, "assemblies", "genomes")
checkm_dir = os.path.join(outDir, "checkM")
checkm_signature = hashlib.sha256(json.dumps({s: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    for s, p in sorted(current_genomes.items()) if not is_lambda_control(expected_organisms[s])}, sort_keys=True).encode()).hexdigest()
checkm_dir = os.path.join(checkm_dir, 'illumina-' + checkm_signature[:16])
# Stage only assemblies produced/verified against this execution's selected reads.
checkm_out, checkm_ver = {}, None
if any(not is_lambda_control(expected_organisms[s]) for s in current_genomes):
    with tempfile.TemporaryDirectory(prefix='illumina-checkm-') as checkm_inputs:
        for sample, genome in current_genomes.items():
            if not is_lambda_control(expected_organisms[sample]):
                (Path(checkm_inputs) / f'{sample}.fasta').symlink_to(Path(genome).resolve())
        checkm_out, checkm_ver = process_data.checkM_stats(
            genomes_dir=checkm_inputs, cpus=cpus, outdir=checkm_dir, logfile=log,
        )

if species_identification_enabled:
    try:
        reconciled_species = reconcile_batch_identifications(
            species_output_root=os.path.join(outDir, "species_identification"),
            database_root=args.species_database_root,
            expected_organisms=expected_organisms,
            platform="Illumina",
            assemblies={
                sample_id: os.path.join(
                    outDir, "assemblies", "genomes", f"{sample_id}.fasta"
                )
                for sample_id in expected_organisms
            },
            gtdb_metadata_tsv=args.species_gtdb_metadata,
            gtdbtk_command_prefix=(
                "conda", "run", "-n", args.species_gtdbtk_conda_env, "gtdbtk"
            ),
            placement_threads=max(1, int(cpus)),
            checkm_contamination_by_sample={
                sample_id: values[1] for sample_id, values in checkm_out.items()
            },
            maximum_checkm_contamination=args.max_contamination,
        )
        update_summary_identifications(temp_qc_summary, reconciled_species)
    except Exception as e:
        time_print(f"Run-wide species-signal review failed: {e}", "Fail")
        logger(log, f"Run-wide species-signal review failed: {e}")
        sys.exit(1)

# Update the summary file with CheckM results
# qc_summary = os.path.join(qc_out, "qc_summary.tsv")
qc_summary = os.path.join(outDir, f"{run_name}_qc_summary.tsv")
time_print('Final summary', "Header")
emit('assembly_qc', 'completed')
logger(log, 'Final summary', "Header")

tool_entries = [
    ("BactiPipe", bactipipe_version),
    ("BactiPipe database", database_identity_label()),
    (args.assembler, current_assembly_version),
]
if not species_identification_enabled:
    tool_entries.append(("KmerFinder", find_organism.kmerfinder_version))
if checkm_ver:
    tool_entries.append(("CheckM", checkm_ver))
tools = OrderedDict(tool_entries)
invoked_fastp = sorted({s['fastp_version'] for s in read_selections.values() if s['fastp_version']})
for sample in expected_organisms:
    invocation = Path(qc_out) / sample / 'illumina_filter' / 'invocation.json'
    if invocation.is_file():
        invoked_fastp.append(json.loads(invocation.read_text())['version'])
if invoked_fastp:
    tools['Fastp'] = ', '.join(sorted(set(invoked_fastp)))
if any(is_lambda_control(organism) for organism in expected_organisms.values()):
    tools["BLASTn"] = command_version(("blastn", "-version"))
    tools["Minimum Lambda identity"] = f"{args.lambda_min_identity:g}%"
    tools["Minimum Lambda genome coverage"] = f"{args.lambda_min_coverage:g}%"
if species_identification_enabled:
    tools["Species identification database"] = identification_database_version(
        args.species_database_root
    )
    tools["skani"] = command_version(("skani",))
    if identification_database_kind(args.species_database_root) == "speciesfinder_kma":
        tools["SpeciesFinder db"] = os.path.basename(database_identity(args.species_database_root).get("database_prefix", "Not recorded"))
        tools["SpeciesFinder"] = python_package_version("SpeciesFinder")
        tools["KMA"] = command_version(("kma",))
    else:
        tools["Kraken2"] = command_version(("kraken2",))
    if args.species_gtdb_metadata:
        tools["GTDB-Tk"] = command_version(
            ("conda", "run", "-n", args.species_gtdbtk_conda_env, "gtdbtk")
        )

tool_lines = ["\n"] # Blank line before tools section
logger(log, database_log_record(args.species_database_root))
for tool, version in tools.items():
    tool_lines.append(f"{tool}>>{version}")
tool_lines.extend([
    "Minimum depth>>30X (Q ≥30); 40X (29 ≤Q <30); 50X (Q <29); identified E. coli: +10X",
    f"Required Mean quality>>{args.minqual:g}",
    f"Read filtering trigger>>Raw mean Q < {args.minqual+1:g}",
    f"Completeness cutoff>>{args.min_completeness:g}%",
    f"Contamination cutoff>>{args.max_contamination:g}%"
])

sep_len = max(len(line) for line in tool_lines)
sep = "=" * sep_len
tools_block = [sep] + tool_lines + [sep]

header += tools_block


process_data.make_summary(
    qc_summary=qc_summary,
    temp_qc_summary=temp_qc_summary,
    header=header,
    checkm_out=checkm_out,
    logfile=log,
    min_completeness=args.min_completeness,
    max_contamination=args.max_contamination,
    defer_pdf=args.defer_qc_pdf,
)

# Selected FASTQs are reproducible intermediates. Keep compact filter reports,
# full-precision metrics and assembly input identities, not duplicate read data.
for sample, selection in read_selections.items():
    if selection['filter_applied']:
        for selected_path in selection['selected_paths']:
            if os.path.isfile(selected_path):
                os.unlink(selected_path)

# Remove temp_dir
if os.path.exists(temp_dir):
    shutil.rmtree(temp_dir)
    time_print(f"Temporary directory removed.")
    logger(log, f"Temporary directory removed.")

# End the timer
end_time = time.time()

# Calculate elapsed time
elapsed_time = end_time - start_time
hours, remainder = divmod(elapsed_time, 3600)
minutes, seconds = divmod(remainder, 60)

# Print the execution time
end_message = f"PIPELINE COMPLETED: Processed {sample_number} samples in {int(hours)} hours, {int(minutes)} minutes, {seconds:.2f} seconds"
time_print(end_message, "Header")
logger(log, end_message, "Header")
