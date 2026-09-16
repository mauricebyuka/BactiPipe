import os
import shutil
import subprocess
from bactipipe.scripts.utils import fastq_metrics, filtlong_with_metrics, merge_gz_fastqs, merge_gz_fastqs_s3, logger as file_logger
from bactipipe.nanopore_quality import filtlong_min_mean_identity

def filtlong_version():
    vercmd = ["filtlong", "--version"]
    try:
        process = subprocess.Popen(vercmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        filtlong_ver = out.decode('utf-8', errors='replace').strip().split()[-1]
    except Exception as e:
        filtlong_ver = "Unknown"
    return filtlong_ver

def qc_nano(
    fastq_file=None,
    raw_folder=None,
    genome_size=None,
    min_avg_quality=15,
    min_read_length=1000,
    min_coverage=100,
    desired_coverage=300,
    output_fastq="trimmed_reads.fastq",
    output_dir=".",
    cpus="8",
    s3=False,
    bucket=None,
    logfile=None
):
    s_name = os.path.basename(output_fastq).replace(".gz", "").replace(".fastq", "")
    # Filtlong expects mean percent correctness, not a Phred-scale number.
    filtlong_quality = str(filtlong_min_mean_identity(min_avg_quality))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    worker_log = os.path.join(output_dir, f"{s_name}.log")
    def log(msg, level="Norm", mode="timestamp"):
        file_logger(worker_log, msg, message_type=level, mode=mode)

    if genome_size is None:
        log(
            f"[{s_name}] Starting size-independent QC; the organism genome size "
            "is not known before assembly."
        )
    else:
        log(f"[{s_name}] Starting QC for {s_name} with genome size {genome_size} and min coverage {min_coverage}X.")

    if s3 and not bucket:
        log("S3 mode requires bucket to be set.")
        return

    try: 
        clutter = []
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_fastq = os.path.join(output_dir, output_fastq)
        final_qual_file = os.path.join(output_dir, "quality_metrics.txt")
        if not os.path.exists(final_qual_file):
            if raw_folder:
                fastq_file = os.path.join(output_dir, os.path.basename(raw_folder) + ".fastq.gz")
                log(f"Processing: {os.path.basename(raw_folder)}--->{os.path.basename(output_fastq)}\n")
                clutter.append(fastq_file)
                log(f"[{s_name}] Consolidating nanopore raw reads into one file...")
                if not os.path.exists(fastq_file):
                    if s3:
                        s3_prefix = f"s3://{bucket}/{raw_folder.strip('/')}/"
                        if not merge_gz_fastqs_s3(s3_prefix, fastq_file):
                            return

                    elif not merge_gz_fastqs(raw_folder, fastq_file):
                        return
                    
                else:
                    log(f"File {fastq_file} already exists. Skipping consolidation.")
            else:
                log(f"Processing: {os.path.basename(fastq_file)} ---> {os.path.basename(output_fastq)}\n")

            # Step 1: Compute raw FASTQ metrics
            log(f"[{s_name}] Computing quality metrics for raw reads...")
            
            try:
                m1 = fastq_metrics(fastq_file, use_external_decompress=True, pigz_threads=int(cpus))
            except Exception as e:
                log(f"[{s_name}] Error computing FASTQ metrics: {e}", "Fail")
                return

            with open(final_qual_file, "w") as f:
                f.write(f"Raw Reads: {m1['reads']:,}\n")
                f.write(f"Raw Total Bases: {m1['total_bases']:,}\n")
                f.write(f"Raw Mean Read Quality: {m1['mean_q_read_np']:.2f}\n")
                f.write(f"Raw Mean Base Quality: {m1['mean_q_base_np']:.2f}\n")
                f.write(f"Raw Mean Quality per Read (Arithmetic): {m1['mean_q_per_read_arith']:.2f}\n")
                f.write(f"Raw Mean Quality per Base (Arithmetic): {m1['mean_q_per_base_arith']:.2f}\n")

            total_bases = m1['total_bases']
            log(f"\t[{s_name}] ---> Raw Total Bases: {total_bases:,}", mode='simple')
            mean_quality = m1['mean_q_read_np']
            log(f"\t[{s_name}] ---> Raw Mean Quality: {mean_quality:.2f}", mode='simple')

            # Trim reads using Filtlong. Unknown organisms intentionally omit
            # size-dependent target-bases/downsampling arguments.
            trimmed_output = output_fastq 

            log(f"[{s_name}] Trimming reads using filtlong (version: {filtlong_version()})...")
            log(f"[{s_name}] Nanopore mean-read Q >= {min_avg_quality}; Filtlong mean identity >= {filtlong_quality}%.")
            if genome_size is None:
                coverage_target = None
                filtlong_command = [
                    "filtlong",
                    "--min_mean_q",
                    filtlong_quality,
                    "--min_length",
                    str(min_read_length),
                    fastq_file,
                ]
            elif total_bases < min_coverage * genome_size:
                min_required_bases = min_coverage * genome_size
                log(f"[{s_name}] Insufficient data for {min_coverage}X coverage. Total bases available: {total_bases:,}. Required: {min_required_bases:,}.")
                log(f"[{s_name}] Keeping the best 90% of reads...")
                filtlong_command = ["filtlong", "--min_mean_q", filtlong_quality, "--keep_percent", "90", "--min_length", str(min_read_length), fastq_file]
                coverage_target = 1
            
            else:
                log(f"[{s_name}] QC: Min quality = {min_avg_quality}, target coverage: {desired_coverage}X coverage")
                coverage_target = desired_coverage
                required_bases = int(coverage_target * genome_size)
                filtlong_command = ["filtlong", "--min_mean_q", filtlong_quality, "--target_bases", str(required_bases), "--min_length", str(min_read_length), fastq_file]

            try:
                os.makedirs(os.path.dirname(trimmed_output), exist_ok=True)
                log(f"[{s_name}] Running filtlong and recalculating metrics...")

                m2 = filtlong_with_metrics(
                    filtlong_cmd=filtlong_command,
                    out_fastq_path=trimmed_output,
                )
            except Exception as e:
                log(f"[{s_name}] Error running filtlong: {e}", "Fail")
                return

            new_mean_quality = m2['mean_q_read_np']
            trim_bases = m2['total_bases']
            log(f"\t[{s_name}] ---> Trimmed Total Bases: {trim_bases:,}", mode='simple')
            log(f"\t[{s_name}] ---> Trimmed Mean Quality: {new_mean_quality:.2f}", mode='simple')

            with open(final_qual_file, "a") as f:
                f.write(f"Trimmed Reads: {m2['reads']:,}\n")
                f.write(f"Trimmed Total Bases: {m2['total_bases']:,}\n")
                f.write(f"Trimmed Mean Read Quality: {m2['mean_q_read_np']:.2f}\n")
                f.write(f"Trimmed Mean Base Quality: {m2['mean_q_base_np']:.2f}\n")
                f.write(f"Trimmed Mean Quality per Read (Arithmetic): {m2['mean_q_per_read_arith']:.2f}\n")
                f.write(f"Trimmed Mean Quality per Base (Arithmetic): {m2['mean_q_per_base_arith']:.2f}\n")

            trim_coverage = int(trim_bases / genome_size) if genome_size else None

            # Check if the new quality meets the threshold
            if new_mean_quality >= min_avg_quality:
                log(f"[{s_name}] Quality threshold met. Average quality: {new_mean_quality:.2f}. Saving data...\n", "Pass")
                # os.rename(f"{trimmed_output}.gz", output_fastq)
            else:
                log(f"[{s_name}] Quality threshold not met. Average quality: {new_mean_quality:.2f}.", "Warn")

            if trim_coverage is None:
                log(
                    f"[{s_name}] Trimmed reads saved to {output_fastq}; final depth will use "
                    f"the assembled genome length. Average quality: {new_mean_quality:.2f}."
                )
            else:
                log(f"[{s_name}] Trimmed reads saved to {output_fastq} with {trim_coverage}X coverage and average quality {new_mean_quality:.2f}.")
            # Cleanup
            log(f"[{s_name}] Cleaning up intermediate files...")
            for item in clutter:
                if os.path.exists(item):
                    if os.path.isdir(item):
                        shutil.rmtree(item)
                    else:
                        os.remove(item)
            if coverage_target is not None and coverage_target < min_coverage:
                log(f"[{s_name}] Unable to achieve {min_avg_quality} average quality with minimum coverage of {min_coverage}X.")
        else:
            log(f"[{s_name}] Quality metrics file already exists. Skipping quality assessment.")

    except Exception as e:
        log(f"[{s_name}] Unexpected error in qc_nano: {e}", "Fail")
        return
