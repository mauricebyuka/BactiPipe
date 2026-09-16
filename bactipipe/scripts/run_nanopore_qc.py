import sys
import os
import csv
import time
import psutil
import shutil
import socket
import argparse
import subprocess
from datetime import datetime
from pathlib import Path
from bactipipe.scripts import nano_qc
from bactipipe.scripts import process_data
from bactipipe.scripts import find_organism
from bactipipe.scripts import sample_failures as sample_failure_audit
from bactipipe.scripts.sample_update import run_qc_sample_update
from bactipipe.scripts.utils import compress_qc_fastqs, logger, pipeheader, excel_reader, time_print, simple_print
from bactipipe.genome_sizes import GenomeSizeCatalog
from bactipipe.assembly_depth import assembly_depth
from bactipipe.lambda_control import is_lambda_control, control_summary_row, assess_control, add_control_arguments
from bactipipe.adaptive_assembly import adaptive_flye
from bactipipe.illumina_parallel import cpu_budget, available_memory_gb, atomic_json
from bactipipe.nanopore_parallel import plan_flye_resources, run_assembly_jobs
from bactipipe.progress import emit
from dataclasses import asdict
from bactipipe.scripts.qc_plan import (
    NanoporeAssemblyJob,
    NanoporeQCJob,
    bounded_worker_count,
    nanopore_qc_failure_reason,
    write_execution_plan,
)
from tqdm.contrib.concurrent import process_map
from collections import OrderedDict
from bactipipe.__version__ import __version__
from bactipipe.reproducibility import database_identity, database_identity_label, database_log_record
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


class CustomHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs['width'] = 100  # Set the width of the help message
        super().__init__(*args, **kwargs)

    def _get_help_string(self, action):
        if action.default is not None and action.default != argparse.SUPPRESS:
            return super()._get_help_string(action)
        return action.help

parser = argparse.ArgumentParser(
    description="Perform quality control on nanopore reads, assemble genomes for reads that pass QC, and check quality for assembled genomes.",
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

optional_args.add_argument(
    "--input-run-name",
    help=(
        "Optional physical source-run directory name. Outputs continue to use "
        "--run_name; intended for controlled reanalysis without copying FASTQs."
    ),
)

required_args.add_argument("-l", "--sample_sheet",
                    help="A two- or three-column file with sample, organism species, and optional barcode. Without a barcode, FASTQs are read from a folder named for the sample.", required=True)

required_args.add_argument('-s', '--source', help='The source directory.', required=True)

optional_args.add_argument("-o", "--outdir",
                    help="Path to output directory. Default: current directory.")
optional_args.add_argument('-t', '--threads', help='Number of CPUs per sample. Default: all available.', type=int)


optional_args.add_argument("--cpus_per_sample", help="Number of CPUs per sample. Must be lower than total allocated CPUs. If not provided,  it will be calculated based on available resources.", type=int)

optional_args.add_argument("-m", "--mincov", help="Minimum genome coverage depth.", default=50, type=int)
optional_args.add_argument("-q", "--minqual", help="Average read quality threshold.", default=15, type=int)
optional_args.add_argument("--min-read-length", help="Minimum read length retained by Filtlong.", default=1000, type=int)
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
optional_args.add_argument ("-a", "--assembler", help="Genome assembler to use.",  choices=["Flye", "Unicycler"], default="Flye")

optional_args.add_argument("--use-s3", help="Use S3 for input. Default: False", action="store_true")

optional_args.add_argument("--bucket", help="S3 bucket name. Required if --use-s3 is set.", type=str, default=None)

optional_args.add_argument("--clean",
                            help="This will delete all fastq samples generated by QC step.",
                            action="store_true")

optional_args.add_argument("-h", "--help",
                           action="help",
                           help="Show this help message and exit.")

add_control_arguments(optional_args)
args = parser.parse_args()

species_identification_enabled = bool(args.species_database_root)
if args.species_gtdb_metadata and not species_identification_enabled:
    parser.error("--species-gtdb-metadata requires --species-database-root")

# Version of BactiPipe
bactipipe_version = f"v{__version__}"

    # Function to run quality control for a single sample
def _process_sample(job):
    sample = job.sample
    output_dir = job.output_dir
    failure_path = os.path.join(output_dir, "analysis_failure.txt")
    if os.path.isfile(failure_path):
        os.remove(failure_path)

    # Run QC with controlled CPU allocation
    if not sample.startswith("#"):
        nano_qc.qc_nano(
            raw_folder=job.raw_folder,
            genome_size=job.genome_size_bp,
            min_avg_quality=job.min_avg_quality,
            min_read_length=job.min_read_length,
            min_coverage=job.min_coverage,
            desired_coverage=job.desired_coverage,
            output_fastq=job.output_fastq,
            output_dir=output_dir,
            cpus=str(job.cpus),
            s3=job.use_s3,
            bucket=job.bucket,
            logfile=os.path.join(output_dir, f"{sample}_qc_nano.log")
        )
        metrics_path = os.path.join(output_dir, "quality_metrics.txt")
        metrics_text = ""
        if os.path.isfile(metrics_path):
            with open(metrics_path, "r", encoding="utf-8", errors="replace") as handle:
                metrics_text = handle.read()
        failure_reason = nanopore_qc_failure_reason(
            metrics_text,
            genome_size_bp=job.genome_size_bp,
            minimum_quality=job.min_avg_quality,
            minimum_coverage=job.min_coverage,
        )
        if not os.path.isfile(os.path.join(output_dir, job.output_fastq)):
            failure_reason = "Read quality filtering did not produce a filtered FASTQ."
        if failure_reason:
            with open(
                failure_path,
                "w",
                encoding="utf-8",
            ) as handle:
                handle.write(failure_reason + "\n")


def process_sample(job):
    """Keep an unexpected worker exception isolated to its sample."""

    sample = job.sample
    emit('read_qc', sample=sample)
    try:
        return _process_sample(job)
    except Exception as exc:
        output_dir = job.output_dir
        os.makedirs(output_dir, exist_ok=True)
        reason = f"Unexpected QC worker failure: {exc}"
        with open(
            os.path.join(output_dir, "analysis_failure.txt"),
            "w",
            encoding="utf-8",
        ) as handle:
            handle.write(reason + "\n")
        return None
    finally:
        emit('read_qc', 'failed' if os.path.isfile(os.path.join(job.output_dir, 'analysis_failure.txt')) else 'completed', sample=sample)

# Function for genome assembly for one sample
if args.assembler:
    assembler = args.assembler.strip().lower()
else:
    assembler = "flye"
def assemble_sample(job):
    sample = job.sample
    # Assemble the genomes of the samples that passed the quality control
    genome = os.path.join(job.assembly_dir, "genomes", f"{sample}.fasta")
    if os.path.exists(genome):
        emit('assembly', 'skipped', sample=sample, metrics={'reused': True})
        return
    emit('assembly', sample=sample, metrics={'threads': job.cpus})
    try:
        if job.assembler == "flye" and job.genome_size_bp is None:
            adaptive_flye(
                reads=job.fastq, output=genome,
                audit_dir=os.path.join(os.path.dirname(job.assembly_dir), "qc_out", sample, "adaptive_assembly"),
                cpus=job.cpus, min_quality=job.min_quality, min_length=job.min_length,
                report=lambda message: logger(job.logfile, f"{sample}: {message}"),
            )
            return
        process_data.assemble(
            sample=sample,
            reads=[job.fastq],
            assembly_dir=job.assembly_dir,
            assembler=job.assembler,
            sequencer="nanopore",
            cpus=job.cpus,
            logfile=job.logfile,
            gsize=job.genome_size_bp,
            single=False,
        )
    except Exception as exc:
        logger(job.logfile, f"Unexpected assembly worker failure for {sample}: {exc}")
        failure_dir = os.path.join(os.path.dirname(job.assembly_dir), "qc_out", sample)
        os.makedirs(failure_dir, exist_ok=True)
        with open(os.path.join(failure_dir, "assembly_failure.txt"), "w") as handle:
            handle.write(str(exc) + "\n")
    return

# Run the pipeline
tech_name = args.name
run_name = args.run_name
source = args.source
source = os.path.join(source, args.input_run_name or run_name)
mincov = args.mincov
minqual = args.minqual
if args.outdir:
    outDir = os.path.join(args.outdir, run_name)
else:
    outDir = os.path.join(os.getcwd(), run_name)


def identify_sample_with_ani(sample, expected_organism, genome, reads, threads):
    """Run the opt-in ANI workflow with candidate-only Kraken evidence."""
    result = identify_assembly(
        sample_id=sample,
        assembly=genome,
        database_root=args.species_database_root,
        output_directory=os.path.join(outDir, "species_identification", sample),
        expected_organism=expected_organism,
        reads=reads,
        platform="Nanopore",
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
            outdir_parent=os.path.dirname(os.path.abspath(outDir)),
            run_name=run_name,
            sample_names=args.only_sample,
            species_database_root=args.species_database_root,
            species_platform="Nanopore",
            species_gtdb_metadata=args.species_gtdb_metadata,
            species_gtdbtk_conda_env=args.species_gtdbtk_conda_env,
            species_placement_threads=max(1, int(args.cpus_per_sample or args.threads or 1)),
        )
    except (OSError, RuntimeError, ValueError) as exc:
        sys.stderr.write(f"BactiPipe sample update failed: {exc}\n")
        sys.exit(1)
    sys.exit(0)


sample_list = os.path.abspath(args.sample_sheet)

log_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
log = os.path.join(outDir, f"{log_time}~{run_name}_qc.log")

date = datetime.now().strftime("%Y-%m-%d")
hostname = socket.gethostname()
def new_func(s):
    s.connect(('10.254.254.254', 1))

try:
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        # doesn't even have to be reachable
        new_func(s)
        ip_address = s.getsockname()[0]
    except Exception as e:
        logger(log, f"Unable to resolve hostname: {e}")
        time_print(f"Unable to resolve hostname: {e}", "Fail")
        ip_address = "Not detected"
    finally:
        s.close()
except Exception as e:
    ip_address = f"Unable to resolve hostname: {e}"
    ip_address = "Unable to resolve hostname"

# Manage computational resources
cpus = cpu_budget(args.threads)


opt_cpus = max(1, cpus // 8)  # type: ignore # Optimal number of CPUs per sample
if opt_cpus == 1:
    defaut_max_samps = 1 # We can only run one sample at a time
else:
    defaut_max_samps = cpus // opt_cpus  # Maximum number of samples to process in parallel

# Calculate available memory
mem = psutil.virtual_memory()
a_mem = mem.available / 1024**3  # available memory in GB
t_mem = mem.total / 1024**3  # total memory in GB
 
allowed_samps = max(1, int(a_mem // 8))  # 8 GB per sample

if allowed_samps < 1:
    logger(log, "Insufficient memory available to run the pipeline.")
    time_print("Insufficient memory available to run the pipeline.", "Fail")
    sys.exit(1)

if allowed_samps < defaut_max_samps:
    pool_size = int(allowed_samps)
else:
    pool_size = int(defaut_max_samps)

if args.cpus_per_sample:
    cpus_per_sample = min(cpus, max(1, int(args.cpus_per_sample)))
else:
    cpus_per_sample = max(1, cpus // pool_size) # type: ignore

# Ensure that the source directory exists
if not os.path.exists(source):
    # print(f"Source directory {source} does not exist.")
    logger(log, f"Source directory {source} does not exist.")
    time_print(f"Source directory {source} does not exist.", "Fail")
    sys.exit(1)

fastq_pass = None
for root, dirs, files in os.walk(source):
    for dir_name in dirs:
        if dir_name == 'fastq_pass':
            fastq_pass = os.path.join(root, 'fastq_pass')
            break
    if fastq_pass:
        break

if fastq_pass:
    source = fastq_pass
else:
    logger(log, f"fastq_pass sub-folder not found in {source}.")
    time_print(f"fastq_pass sub-folder not found in {source}.", "Fail")
    sys.exit(1)

header_info = pipeheader("Oxford Nanopore", date, tech_name, hostname, ip_address, run_name, sample_list, source, outDir, cpus)

header = header_info[0]
run_info = header_info[1]
for line in header:
    simple_print(line)
    logger(log, line, mode="simple")
for line in run_info:
    simple_print(line)
    logger(log, line, mode="simple")

qc_out = os.path.join(outDir, "qc_out")
if not os.path.exists(qc_out):
    os.makedirs(qc_out)


# Ensure that the sample list if formatted correctly
if not os.path.exists(sample_list):
    # print(f"Sample list file {sample_list} does not exist.")
    logger(log, f"Sample list file {sample_list} does not exist.")
    time_print(f"Sample list file {sample_list} does not exist.", "Fail")
    sys.exit(1)
else:
    if sample_list.endswith('.xlsx') or sample_list.endswith('.xls'):
        # If the sample list is an Excel file, read it
        try:
            sample_info = excel_reader(sample_list)
        except Exception as e:
            logger(log, f"Error reading Excel file {sample_list}: {e}")
            time_print(f"Error reading Excel file {sample_list}: {e}", "Fail")
            sys.exit(1)
    else:
        with open(sample_list, 'r') as sampL:
            sample_info = sampL.readlines()

        first_line = sample_info[0].strip().split('\t')
        if len(first_line) not in {2, 3}:
            logger(log, "Sample list file must have sample, organism, and an optional barcode column.")
            time_print("Sample list file must have sample, organism, and an optional barcode column.", "Fail")
            sys.exit(1)
        if len(first_line) == 2:
            sample_info = [
                f"{line.strip()}\t{line.strip().split(chr(9))[0]}\n"
                for line in sample_info
                if line.strip()
            ]

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
        logger(log, message)
        time_print(message, "Fail")
        sys.exit(1)
    sample_info = selected_lines

# Resolve the controlled organism and its size exactly once per sample.  The
# same immutable resolution is used by QC, assembly, and final reporting.
genome_size_catalog = GenomeSizeCatalog.from_package()

# Validate the sample sheet
bad_organisms = []
bad_samples = []
bad_barcodes = []
sample_resolutions = {}
sample_failures = {}
expected_organisms = {}
source_entries = set(os.listdir(source))


# with open(sample_list, 'r') as sampL:
for line in sample_info:
    if not line.strip() or line.startswith("#"):
        continue
    sample, organism, barcode = line.strip().split('\t')
    if sample.strip().casefold() in {"sample", "sample_id"}:
        continue
    expected_organisms[sample] = organism
    resolution = genome_size_catalog.resolve(organism)
    sample_resolutions[sample] = resolution
    if resolution.status in {"unsupported", "missing"}:
        reason = (
            f"Organism [{organism.strip()}] is not in the controlled BactiPipe "
            "organism list or was left blank. Continuing in Unknown mode."
        )
        logger(log, f"WARNING: Sample {sample}: {reason}")
        time_print(f"WARNING: Sample {sample}: {reason}", "Fail")
        bad_organisms.append(f"{sample}:{organism}")
    
    if sample not in source_entries and barcode.replace("NB", "barcode") not in source_entries:
        logger(log, f"WARNING: Sample {sample} does not have corresponding fastq files.")
        time_print(f"WARNING: Sample {sample} does not have corresponding fastq files.", "Fail")
        bad_samples.append(sample)
        sample_failure_audit.record_failure(
            sample_failures,
            sample,
            "input",
            "No corresponding FASTQ files were found.",
        )
    if barcode == sample:
        pass
    elif "NB" not in barcode and "barcode" not in barcode:
        logger(log, f"Barcode for sample {sample} is not in the correct format.")
        time_print(f"Barcode for sample {sample} is not in the correct format.", "Fail")
        bad_barcodes.append(f"{sample}:{barcode}")
    elif "NB" in barcode and len(barcode) != 4:
        logger(log, f"Barcode for sample {sample} is not in the correct format.")
        time_print(f"Barcode for sample {sample} is not in the correct format.", "Fail")
        bad_barcodes.append(f"{sample}:{barcode}")
    elif "barcode" in barcode and len(barcode) != 9:
        logger(log, f"Barcode for sample {sample} is not in the correct format.")
        time_print(f"Barcode for sample {sample} is not in the correct format.", "Fail")
        bad_barcodes.append(f"{sample}:{barcode}")

if bad_barcodes:
    logger(log, "Some sample-sheet barcodes are invalid; affected samples will be skipped:")
    time_print("Some sample-sheet barcodes are invalid; affected samples will be skipped:", "Fail")
    logger(log, f" --> Invalid barcodes: {', '.join(bad_barcodes)}")
    time_print(f" --> Invalid barcodes: {', '.join(bad_barcodes)}")
    for item in bad_barcodes:
        sample = item.split(":", 1)[0]
        if sample not in bad_samples:
            bad_samples.append(sample)
        sample_failure_audit.record_failure(
            sample_failures,
            sample,
            "input",
            "The barcode value is not in a supported format.",
        )

if bad_organisms:
    time_print("Unsupported or blank organisms will use Unknown mode.")
    logger(log, "Unsupported or blank organisms will use Unknown mode.")

if bad_samples:
    b_s = ", ".join(bad_samples)
    time_print(f"WARNING: These samples don't have corresponding fastq data. They will be skipped: {b_s}.", "Fail")
    logger(log, f"WARNING: These samples don't have corresponding fastq data. They will be skipped: {b_s}.")

sample_lines = []
for line in sample_info:
    if not line.strip() or line.startswith("#"):
        continue
    sample = line.strip().split('\t')[0]
    if sample.strip().casefold() not in {"sample", "sample_id"} and sample not in bad_samples:
        sample_lines.append(line.strip())

sample_number = len(sample_lines)

if sample_number == 0:
    time_print("No valid samples to process. Exiting pipeline.", "Fail")
    logger(log, "No valid samples to process. Exiting pipeline.")
    sys.exit(1)

pool_size = bounded_worker_count(
    sample_count=sample_number,
    total_cpus=cpus,
    cpus_per_sample=cpus_per_sample,
    memory_worker_limit=allowed_samps,
)

write_execution_plan(
    os.path.join(outDir, "sample_execution_plan.tsv"),
    [
        (sample, expected_organisms[sample], sample_resolutions[sample])
        for sample in expected_organisms
    ],
    platform="Nanopore",
    min_quality=minqual,
    min_read_length=args.min_read_length,
    min_coverage=mincov,
    desired_coverage=300,
    catalog_checksum=genome_size_catalog.checksum,
)

qc_jobs = []
for line in sample_lines:
    sample, organism, barcode = line.split("\t")
    normalized_barcode = barcode.replace("NB", "barcode")
    resolution = sample_resolutions[sample]
    qc_jobs.append(
        NanoporeQCJob(
            sample=sample,
            organism_entered=organism,
            organism_canonical=resolution.canonical_name or "Unknown",
            barcode=normalized_barcode,
            raw_folder=os.path.join(source, normalized_barcode),
            output_dir=os.path.join(qc_out, sample),
            output_fastq=f"{sample}.fastq",
            genome_size_bp=resolution.genome_size_bp,
            min_avg_quality=minqual,
            min_read_length=args.min_read_length,
            min_coverage=mincov,
            desired_coverage=300,
            cpus=cpus_per_sample,
            use_s3=args.use_s3,
            bucket=args.bucket if args.use_s3 else None,
        )
    )

# Compute pool size (number of parallel samples)
print(f"\nRunig QC on samples in run {run_name}.\n")
logger(log, f"Running QC on samples in run {run_name}.\n")
print(f"  --> Total CPUs available: {cpus}")
logger(log, f"  --> Total CPUs available: {cpus}")
print(f"  --> Total memory available: {int(a_mem)} GB out of {int(t_mem)} GB.")
logger(log, f"  --> Total memory available: {int(a_mem)} GB out of {int(t_mem)} GB.")
print(f"  --> Using {pool_size} parallel samples, each using {cpus_per_sample} CPU(s).\n")
logger(log, f"  --> Using {pool_size} parallel samples, each using {cpus_per_sample} CPU(s).")


start_time = time.time()
logger(log, "PIPELINE STARTED", "Header")
time_print("PIPELINE STARTED", "Header")


qc_start_msg = "Performing Quality Control"
emit('inputs', 'completed')
emit('read_qc', total=len(qc_jobs))
time_print(qc_start_msg, "Header")
logger(log, qc_start_msg, "Header")

time_print(f"\nTotal number of samples to be processed: {sample_number}\n")
print()
logger(log, f"Total number of samples to be processed: {sample_number}\n")


# QC stage (multiprocessing) ----
q_bar_fmt = '{l_bar} {bar:40} {n_fmt}/{total_fmt} [{elapsed}{postfix}]'
process_map(
    process_sample,
    qc_jobs,
    max_workers=pool_size,
    bar_format=q_bar_fmt,
    chunksize=1,
    desc="---> Processing fastq files: ",
    colour="yellow",
    disable=False,
)

for line in sample_lines:
    sample = line.split('\t', 1)[0]
    failure_path = os.path.join(qc_out, sample, "analysis_failure.txt")
    if os.path.isfile(failure_path):
        with open(failure_path, "r", encoding="utf-8", errors="replace") as handle:
            reason = handle.read().strip()
        sample_failure_audit.record_failure(
            sample_failures, sample, "qc", reason
        )

# Count the number of samples that passed QC
passed_samples = [
    line.split('\t', 1)[0]
    for line in sample_lines
    if not os.path.isfile(
        os.path.join(qc_out, line.split('\t', 1)[0], "analysis_failure.txt")
    )
    and os.path.exists(
        os.path.join(qc_out, line.split('\t', 1)[0], "quality_metrics.txt")
    )
]

passed_number = len(passed_samples)
emit('read_qc', 'completed', total=len(qc_jobs), finished=len(qc_jobs), failed=len(qc_jobs) - passed_number)
print()

if len(passed_samples) == 0:
    time_print("No samples passed QC. Exiting pipeline.", "Fail")
    logger(log, "No samples passed QC. Exiting pipeline.")
    sys.exit(1)

# Print samples that failed QC
failed_samples = [sample for sample in os.listdir(qc_out) if not os.path.exists(os.path.join(qc_out, sample, "quality_metrics.txt")) and sample != "temp_qc_summary.tsv"]

if failed_samples:
    print()
    failures = ".".join(failed_samples)
    fail_msg = f'QC failed for samples: {failures}.'
    time_print(fail_msg, "Fail")
    logger(log, fail_msg)

# Assemble genomes for samples that passed QC
assembly_start_msg = "Genome Assembly for samples that passed QC" 
time_print(assembly_start_msg, "Header")
emit('assembly', total=passed_number)
logger(log, assembly_start_msg, "Header")

fastq_files = []
for sample in sorted(os.listdir(qc_out)):
    sample_dir = os.path.join(qc_out, sample)
    if not os.path.isdir(sample_dir) or sample not in passed_samples:
        continue
    p_fastq = os.path.join(qc_out, sample, f"{sample}.fastq")
    p_gz = os.path.join(qc_out, sample, f"{sample}.fastq.gz")
    if os.path.exists(p_fastq):
        fastq_files.append((sample, p_fastq))
    elif os.path.exists(p_gz):
        fastq_files.append((sample, p_gz))

pass_msg = f"Available for assembly: {len(fastq_files)} of {sample_number} samples."
time_print(pass_msg, "Pass")
logger(log, pass_msg)

genomes_dir = os.path.join(outDir, "assemblies", "genomes")
if not os.path.exists(genomes_dir):
    os.makedirs(genomes_dir)

existing_genomes = [
    f.split(".")[0]
    for f in os.listdir(genomes_dir)
    if f.endswith(".fasta")
    and f.split(".")[0] in passed_samples
    and os.path.getsize(os.path.join(genomes_dir, f)) > 0
]
assembly_input = []
num_to_assemble = 0

if len(existing_genomes) == 0:
    num_to_assemble = len(fastq_files)
    assembly_input = fastq_files

    print(f"  ---> Genomes for {len(fastq_files)} samples will be assembled.")
    logger(log, f"  ---> Genomes for {len(fastq_files)} samples will be assembled.")

elif len(existing_genomes) > 0 and len(existing_genomes) < len(fastq_files):
    num_to_assemble = len(fastq_files) - len(existing_genomes)
    assembly_input = [item for item in fastq_files if item[0] not in existing_genomes]
    time_print(f"Genomes for {len(existing_genomes)} samples have been assembled already.")
    logger(log, f"Genomes for {len(existing_genomes)} samples have been assembled already.")

    simple_print(f"---> Genomes for {num_to_assemble} samples will be assembled.")
    logger(log, f"---> Genomes for {num_to_assemble} samples will be assembled.")

elif len(existing_genomes) == passed_number:
    num_to_assemble = 0
    assembly_input = []
    time_print("Genomes for all samples have been assembled already. Skiping assembly!")
    logger(log, "Genomes for all samples have been assembled already. Skiping assembly!")

if not assembly_input:
    if existing_genomes:
        time_print("Genomes for all available samples have been assembled already. Skipping assembly.")
        logger(log, "Genomes for all available samples have been assembled already. Skipping assembly.")
    else:
        time_print("No samples to assemble. Exiting pipeline.", "Fail")
        logger(log, "No samples to assemble. Exiting pipeline.")
        sys.exit(1)

if assembly_input:
    if assembler == 'flye':
        assembly_resources = plan_flye_resources(num_to_assemble, cpus,
            requested_threads=args.cpus_per_sample, memory_gb=available_memory_gb())
        assembly_pool_size = assembly_resources.jobs
        assembly_cpus = assembly_resources.threads_per_sample
        atomic_json(Path(outDir) / 'nanopore_assembly_resources.json', asdict(assembly_resources))
    else:
        assembly_cpus = cpus_per_sample if args.cpus_per_sample else max(1, cpus // num_to_assemble)
        if assembler == 'unicycler':
            assembly_cpus = min(cpus, max(assembly_cpus, 10))
        assembly_pool_size = bounded_worker_count(sample_count=num_to_assemble, total_cpus=cpus,
            cpus_per_sample=assembly_cpus, memory_worker_limit=allowed_samps)
    assembly_dir = os.path.join(outDir, "assemblies")
    assembly_jobs = [
        NanoporeAssemblyJob(
            sample=sample,
            fastq=fastq,
            genome_size_bp=sample_resolutions[sample].genome_size_bp,
            assembly_dir=assembly_dir,
            assembler=assembler,
            cpus=assembly_cpus,
            logfile=log,
            min_quality=minqual,
            min_length=args.min_read_length,
        )
        for sample, fastq in assembly_input
    ]

    stdout = subprocess.run([assembler, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = (stdout.stdout or stdout.stderr).decode('utf-8', errors='replace')
    assembly_version = out.split()[-1] if out.split() else "Unknown"
    if not assembly_version.lower().startswith("v"):
        assembly_version = "v" + assembly_version

    assembler_message = f"Using {assembler.capitalize()} (version: {assembly_version})."

    simple_print(f"  ---> {assembler_message}")
    logger(log, f"  ---> {assembler_message}")

    simple_print(f"  ---> {assembly_cpus} CPUs will be used per sample.")
    logger(log, f"  ---> {assembly_cpus} CPUs will be used per sample.")
    allocation = (f"Assembly allocation: {assembly_pool_size} concurrent samples x {assembly_cpus} threads "
                  f"= {assembly_pool_size * assembly_cpus}/{cpus} allocated CPUs; {num_to_assemble} samples total.")
    emit('assembly', total=num_to_assemble, metrics={'concurrent_samples': assembly_pool_size,
         'threads_per_sample': assembly_cpus, 'allocated_cpus': cpus})
    time_print(allocation)
    logger(log, allocation)
    sys.stdout.flush()

    def report_assembly_progress(message):
        time_print(message)
        logger(log, message)
        sys.stdout.flush()
    run_assembly_jobs(assembly_jobs, workers=assembly_pool_size, worker=assemble_sample,
        report=report_assembly_progress, state_path=Path(outDir) / 'nanopore_assembly_progress.json')

# Count the number of samples that have been assembled
assembled_samples = [f for f in os.listdir(genomes_dir) if f.endswith(".fasta") and os.path.getsize(os.path.join(genomes_dir, f)) > 0]

assembled_number = len(assembled_samples)
assembly_msg = f"Genome assembly completed for {assembled_number} of {passed_number} samples.\n"
print()
time_print(assembly_msg, "Pass")
emit('assembly', 'completed', finished=passed_number,
     failed=sum(not os.path.isfile(os.path.join(genomes_dir, sample + '.fasta')) for sample in passed_samples))
logger(log, assembly_msg)

# Print samples that failed assembly
failed_assembly = [sample for sample in passed_samples if not os.path.exists(os.path.join(genomes_dir, f"{sample}.fasta")) or os.path.getsize(os.path.join(genomes_dir, f"{sample}.fasta")) == 0]
if failed_assembly:
    failed_assembly_msg = f'Assembly failed for samples: {", ".join(failed_assembly)}.\n'
    time_print(failed_assembly_msg, "Fail")
    logger(log, failed_assembly_msg, "Fail")
    for sample in failed_assembly:
        failure_path = os.path.join(outDir, "qc_out", sample, "assembly_failure.txt")
        assembly_reason = "Genome assembly did not produce a usable FASTA file."
        if os.path.isfile(failure_path):
            with open(failure_path) as handle:
                assembly_reason = handle.read().strip() or assembly_reason
        sample_failure_audit.record_failure(
            sample_failures,
            sample,
            "assembly",
            assembly_reason,
        )

if len(failed_assembly) == passed_number:
    time_print("No samples were assembled. Exiting pipeline.", "Fail")
    logger(log, "No samples were assembled. Exiting pipeline.")
    sys.exit(1)

kmer_message = "Verifying Taxonomic Identity"
time_print(kmer_message, "Header")
emit('identification', total=sum(sample not in bad_samples for sample in expected_organisms))
logger(log, kmer_message, "Header")

# Find the organism for each sample
temp_qc_summary = os.path.join(qc_out, "temp_qc_summary.tsv")
reported_samples = set()
with(open(temp_qc_summary , 'w')) as qc_sum:
    writer = csv.writer(qc_sum, dialect='excel-tab')
    writer.writerow(["Sample",  "Mean quality", "qc_verdict", "Expected organism", "Identified organism", "% Match", "Coverage", "min_cov", "cov_verdict", "tax_confirm"])

    for line in sample_info:
        if line.startswith("#"):
            continue
        sample, organism, barcode = line.strip().split('\t')
        if sample.strip().casefold() in {"sample", "sample_id"}:
            continue
        if sample in bad_samples:
            continue

        emit('identification', sample=sample)
        resolution = sample_resolutions[sample]
        if not resolution.is_known:
            sys_organism = "unknown"
            genome_size = None
        else:
            sys_organism = resolution.canonical_name
            genome_size = resolution.genome_size_bp

        genome = os.path.join(outDir, "assemblies", "genomes", f"{sample}.fasta")
        if not os.path.exists(genome) or os.path.getsize(genome) == 0:
            time_print(f"WARNING: Sample {sample} does not have a valid genome file. Skipping.", "Fail")
            logger(log, f"WARNING: Sample {sample} does not have a valid genome file. Skipping.")
            sample_failure_audit.record_failure(
                sample_failures,
                sample,
                "assembly",
                "Genome assembly did not produce a usable FASTA file.",
            )
            continue
        # Find quality metrics:
        qc_file = os.path.join(qc_out, sample, "quality_metrics.txt")

        with open(qc_file, 'r') as qcF:
            lines = qcF.readlines()
        
        # Extract total bases
        bases_line = next((l for l in lines if "Trimmed Total Bases:" in l), None)
        if bases_line is None:
            bases_line = next((l for l in lines if "Raw Total Bases:" in l), None)
        if bases_line is None:
            logger(log, f"Total bases line not found in QC metrics for sample {sample}")
            time_print(f"Total bases line not found in QC metrics for sample {sample}", "Fail")
            sample_failure_audit.record_failure(
                sample_failures, sample, "qc", "QC metrics are incomplete."
            )
            continue
        try:
            total_bases = int(float(bases_line.split(":")[1].strip().replace(",", "")))
        except (IndexError, TypeError, ValueError) as e:
            sample_failure_audit.record_failure(
                sample_failures, sample, "qc", f"QC total bases could not be parsed: {e}"
            )
            continue

        # Extract mean read quality
        qual_line = next((l for l in lines if "Trimmed Mean Read Quality:" in l), None)
        if qual_line is None:
            qual_line = next((l for l in lines if "Raw Mean Read Quality:" in l), None)
        if qual_line is None:
            logger(log, f"Mean read quality line not found in QC metrics for sample {sample}")
            time_print(f"Mean read quality line not found in QC metrics for sample {sample}", "Fail")
            sample_failure_audit.record_failure(
                sample_failures, sample, "qc", "QC metrics are incomplete."
            )
            continue
        try:
            avqc = float(qual_line.split(":")[1].strip())
        except (IndexError, TypeError, ValueError) as e:
            sample_failure_audit.record_failure(
                sample_failures, sample, "qc", f"Mean read quality could not be parsed: {e}"
            )
            continue

        if avqc < minqual:
            qc_verdict = "Fail"
        else:
            qc_verdict = "Pass"

        coverage, cov_verdict = assembly_depth(
            total_bases, genome, mincov,
            audit_path=os.path.join(qc_out, sample, "depth_metrics.json"),
        )
        cov_display = f"{coverage:.2f}X" if coverage is not None else "N/A"

        if is_lambda_control(organism):
            control_evidence = assess_control(genome, os.path.join(qc_out, sample, "lambda_control"),
                min_identity=args.lambda_min_identity, min_coverage=args.lambda_min_coverage,
                report=lambda message: logger(log, message))
            writer.writerow(control_summary_row(sample, avqc, qc_verdict, organism, cov_display, mincov, cov_verdict, control_evidence))
            if control_evidence["status"] == "Error":
                sample_failure_audit.record_failure(sample_failures, sample, "control_confirmation", control_evidence["reason"])
            reported_samples.add(sample)
            continue

        if species_identification_enabled:
            filtered_reads = [
                path
                for path in (
                    os.path.join(qc_out, sample, f"{sample}.fastq"),
                    os.path.join(qc_out, sample, f"{sample}.fastq.gz"),
                )
                if os.path.isfile(path)
            ][:1]
            try:
                summary_fields = identify_sample_with_ani(
                    sample, sys_organism, genome, filtered_reads, cpus_per_sample
                )
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
            identified_org = summary_fields.coverage_organism
            best_other_org = None
        else:
            try:
                kmer_result = find_organism.find_species_with_kmrf(s_name=sample, lab_species=sys_organism, genome=genome, dataOut=outDir, org_type="bacteria", logfile=log)
                if not kmer_result:
                    sample_failure_audit.record_failure(
                        sample_failures,
                        sample,
                        "organism_identification",
                        "Organism identification did not return a result.",
                    )
                    continue
                hit, tax_confirm, possibilities = kmer_result
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
                    best_percent = hit.split(" (")[1].strip(")")
                    identified_org = best_org
                else:
                    best_other_hit = hit.split(" --- ")[0].split(": ")[1]
                    best_other_org = best_other_hit.split(" (")[0]
                    best_other_percent = best_other_hit.split(" (")[1].strip(")")
                    identified_org = best_other_org
            except (AttributeError, IndexError, TypeError, ValueError) as e:
                sample_failure_audit.record_failure(
                    sample_failures,
                    sample,
                    "organism_identification",
                    f"Organism-identification output could not be parsed: {e}",
                )
                continue

        logger(log, f"Sample {sample}: assembly-based depth {cov_display}; {cov_verdict}")

        if species_identification_enabled or best_org:
            writer.writerow([sample, f'{avqc:.2f}', qc_verdict, organism, best_org, best_percent, cov_display, mincov, cov_verdict, tax_confirm])
        else:
            writer.writerow([sample, f'{avqc:.2f}', qc_verdict, organism, f'Closest: {best_other_org}', best_other_percent, cov_display, mincov, cov_verdict, tax_confirm])
        reported_samples.add(sample)
        emit('identification', 'completed', sample=sample)

sample_failure_audit.publish_failures(
    outDir,
    temp_qc_summary,
    expected_organisms,
    sample_failures,
    required_depth=mincov,
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
checkm_out, checkm_ver = process_data.checkM_stats(
    genomes_dir=genomes_dir, cpus=cpus, outdir=checkm_dir, logfile=log,
    exclude_samples=[sample for sample, organism in expected_organisms.items() if is_lambda_control(organism)],
)

if species_identification_enabled:
    try:
        reconciled_species = reconcile_batch_identifications(
            species_output_root=os.path.join(outDir, "species_identification"),
            database_root=args.species_database_root,
            expected_organisms=expected_organisms,
            platform="Nanopore",
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
            placement_threads=max(1, int(cpus_per_sample)),
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

# Final summary
qc_summary = os.path.join(outDir, f"{run_name}_qc_summary.tsv")
time_print('Final summary', "Header")
emit('assembly_qc', 'completed')
logger(log, 'Final summary', "Header")

tool_entries = [
    ("BactiPipe", bactipipe_version),
    ("BactiPipe database", database_identity_label()),
    ("Filtlong", nano_qc.filtlong_version()),
    (args.assembler, assembly_version),
]
if not species_identification_enabled:
    tool_entries.append(("KmerFinder", find_organism.kmerfinder_version))
tool_entries.append(("CheckM", checkm_ver))
tools = OrderedDict(tool_entries)
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
    f"Minimum depth>>{mincov:g}X",
    f"Required Mean quality>>{args.minqual}",
    f"Minimum read length>>{args.min_read_length} bp",
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

# Delete fastq samples from the QC step
if args.clean:
    for sample in os.listdir(qc_out):
        sample_dir = os.path.join(qc_out, sample)
        if os.path.isdir(sample_dir):
            for file in os.listdir(sample_dir):
                gonome_file = os.path.join(genomes_dir, f"{sample}.fasta")
                if file.endswith(".fastq") and os.path.exists(gonome_file) and os.path.getsize(gonome_file) > 0:
                    os.remove(os.path.join(sample_dir, file)) 
    time_print("Fastq files generated during QC step have been deleted to save space.", "Pass")
    logger(log, "Fastq files generated during QC step have been deleted to save space.")
else:

    def _parent_log_adapter(message: str, level: str = "Norm"):
        logger(log, message, message_type=level)

    # compress uncompressed FASTQs at the end of the pipeline
    summary = compress_qc_fastqs(
        qc_out=qc_out,
        genomes_dir=genomes_dir,
        cpus_per_sample=int(cpus_per_sample),
        require_assembly=False,
        recursive=False,
        remove_source=True,
        logger_fn=_parent_log_adapter,
    )
    logger(log, f"Compression summary: {summary}")

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
print("\n")
