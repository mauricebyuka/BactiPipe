from __future__ import annotations
import os
import re
import sys
import boto3
import shutil
import logging
import textwrap3
import subprocess
import pandas as pd
from tqdm import tqdm
from datetime import datetime
from colorama import Fore, Style, init
from typing import Iterable, Callable, Optional, Dict, List, Tuple, Union
from concurrent.futures import ThreadPoolExecutor, as_completed
import io, gzip, threading, math
from Bio import SeqIO
from bactipipe.__version__ import __version__
from bactipipe.reproducibility import database_identity_label
import csv
import socket
from pathlib import Path

# ReportLab imports
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib import colors
from reportlab.lib.units import inch, cm
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.platypus import SimpleDocTemplate, Paragraph, Table, TableStyle, Spacer
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT

bactipipe_version = f"v{__version__}"

_LEVEL_MAP = {
    "Norm": logging.INFO,
    "Pass": logging.INFO,
    "Header": logging.INFO,
    "Warn": logging.WARNING,
    "Fail": logging.ERROR,
}

def stringwraper(text, width, s_type):
    lines = textwrap3.wrap(text, width=width, break_long_words=True)
    new_lines = []
    if s_type == "plain":
        new_lines = lines
    elif s_type == "command":
        for line in lines[:-1]:
            new_lines.append(line + "\\")
        new_lines.append(lines[-1])
    return new_lines

def time_logger(logfile, message, message_type="", s_type='plain'):
    # logfile is an open file object
    # Adds a timestamp to the message and writes it to the logfile
    tstamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    # message = f"{tstamp} - {message}"
    messag_len = len(tstamp) + 2
    indent = " " * messag_len
    wrap_width = 70#70 - messag_len

    lines = stringwraper(message, width=wrap_width, s_type=s_type)
    if message_type == 'Header':
        logfile.write('\n')
        logfile.write(f'{tstamp} - {message}\n')
        logfile.write("="*len(f'{tstamp} - {message}') + '\n\n')
    elif message_type in ['Fail', 'Pass', 'Norm']:
        logfile.write(f'{tstamp} - {lines[0]}\n')
        for line in lines[1:]:
            logfile.write(indent + line + '\n')
    else:
        logfile.write(message + '\n')

def simple_logger(logfile, message, message_type="", s_type='plain'):
    # This allows for the logging without a timestamp
    if message_type == 'Header':
        logfile.write('\n')
        logfile.write(message + '\n')
        logfile.write("="*len(message) + '\n\n')
    else:
    	logfile.write(message + '\n')

def logger(logfile, message="", message_type="Norm", mode="timestamp", s_type='plain'):
    # Mode options: timestamp, simple
    logfile = os.path.abspath(logfile)
    if not os.path.exists(os.path.dirname(logfile)):
        os.makedirs(os.path.dirname(logfile))
    if mode == "timestamp":
        with open(logfile, 'a') as log:
            time_logger(logfile=log, message=message, message_type=message_type, s_type=s_type)
    elif mode == "simple":
        with open(logfile, 'a') as log:
            simple_logger(logfile=log, message=message, message_type=message_type, s_type=s_type)

def time_print(message, message_type="", s_type='plain'):
    tstamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    # message = f"{tstamp} - {message}"
    messag_len = len(tstamp) + 2
    indent = " " * messag_len
    wrap_width = 140 - messag_len

    lines = stringwraper(message, width=wrap_width, s_type=s_type)
  
    if message_type == 'Fail':
        print(Fore.RED + f'{tstamp} - {lines[0]}' + Style.RESET_ALL)  
        for line in lines[1:]:
            print(Fore.RED + indent + line + Style.RESET_ALL)
        # print('\n')

    elif message_type == 'Pass':
        print(Fore.GREEN + f'{tstamp} - {lines[0]}' + Style.RESET_ALL)  
        for line in lines[1:]:
            print(Fore.GREEN + indent + line + Style.RESET_ALL)
        # print('\n')

    elif message_type == 'Header':
        print('')
        print("\033[1m" + Fore.YELLOW + f'{tstamp} - {message}' + Style.RESET_ALL)  
        print("\033[1m" + Fore.YELLOW + "="*len(f'{tstamp} - {message}') + '\n' + Style.RESET_ALL)  
    else:
        print(f'{tstamp} - {lines[0]}')
        for line in lines[1:]:
            print(indent + line )
            # print('\n')

def simple_print(message, message_type="", s_type='plain'):
    if message_type == 'Fail':
        print(Fore.RED + message.lstrip() + Style.RESET_ALL)  

    elif message_type == 'Pass':
        print(Fore.GREEN + message.lstrip() + Style.RESET_ALL)  

    elif message_type == 'Header':
        print('\n')
        print("\033[1m" + Fore.YELLOW + message.lstrip() + Style.RESET_ALL)  # Green for positive
        print("\033[1m" + Fore.YELLOW + "="*len(message.lstrip()) + '\n' + Style.RESET_ALL)  # Green for positive
    else:
        print(message)  


def check_required(programs, python_packages):
    missing_programs = []
    missing_packages = []
    
    for program in programs:
        result = subprocess.run(['which', program], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            missing_programs.append(program)
    
    for package in python_packages:
        try:
            __import__(package)
        except ImportError:
            missing_packages.append(package)
    
    if missing_programs or missing_packages:
        if missing_programs:
            print("The following required programs are not installed:")
            for program in missing_programs:
                print(f"- {program}")
        if missing_packages:
            print("The following required Python packages are not installed:")
            for package in missing_packages:
                print(f"- {package}")
        print("Please install the missing programs and packages and try again.")
        sys.exit(1)

def pipeheader(platform, date, tech_name, hostname, ip_address, run_name, sample_list, raw_reads, outDir, cpus):
    header_lines = [
        "Bacteria WGS QC Pipeline Results",
        f"Platform: {platform}",
        f"Run name: {run_name}",
        f"Date: {date}",
        f"Analysis by: {tech_name}",
        f"Server Host Name: {hostname}",
    ]
    sep_len = max(len(line) for line in header_lines)
    sep = "=" * sep_len
    header = [sep] + header_lines + [sep]

    run_info = textwrap3.dedent(f'''
    Run Parameters:
    ===============
    \tRun name: {run_name}
    \tSample sheet: {sample_list}
    \tInput directory (fastq reads): {raw_reads}
    \tOutput directory (results): {os.path.abspath(outDir)}
    \tNumber of threads: {cpus}
    ''')
    run_info = run_info.split('\n')
    return header, run_info

def excel_reader(filepath: str, platform: str = "nanopore"):
    # — normalize function treats NaN as empty and collapses whitespace/linebreaks —
    def normalize(x):
        if pd.isna(x):
            return ""
        s = str(x).strip().strip('"')
        return re.sub(r'\s+', ' ', s)

    plat = platform.lower()
    if plat == "nanopore":
        needed = [
            "Accession Number",
            "Host and Source",
            "Bacteria Species",
            "gDNA Prep ID",
        ]
        idx_col = "Native Barcode Index"
        pattern = re.compile(r"^(?:NB\d{2}|barcode\d{2})$", re.IGNORECASE)
        header_lines = 1

    elif plat == "illumina":
        needed = [
            "Well",
            "Accession Number",
            "Host and Source",
            "Bacteria Species",
            "gDNA Prep ID",
            "UD Set A Index Well",
        ]
        idx_col = "Well"
        pattern = re.compile(r"^[A-H](?:0[1-9]|1[0-2])$")
        header_lines = 2

    else:
        raise ValueError("`platform` must be 'nanopore' or 'illumina'")

    # 1) Read every cell as raw (no header)
    raw = pd.read_excel(filepath, sheet_name=0, header=None, engine="openpyxl")
    n_rows, n_cols = raw.shape

    # 2) Find the header block
    header_row = None
    header_vals = None

    for i in range(n_rows - header_lines + 1):
        if header_lines == 1:
            row_vals = [normalize(x) for x in raw.iloc[i]]
            row_vals = ['Native Barcode Index' if x == 'NB Index' else x for x in row_vals]
        else:
            top = [normalize(x) for x in raw.iloc[i]]
            bot = [normalize(x) for x in raw.iloc[i + 1]]
            row_vals = [top[j] or bot[j] for j in range(n_cols)]
        if all(col in row_vals for col in needed):
            header_row = i
            header_vals = row_vals
            if plat == "nanopore" and idx_col in row_vals:
                needed.append(idx_col)
            if "Other Bacteria Species" in row_vals:
                needed.append("Other Bacteria Species")
            break

    if header_row is None:
        raise ValueError(f"Could not locate the header block for platform '{platform}'.")

    # 3) Now safely compute where data starts
    data_start = header_row + header_lines

    # 4) Map each needed name to its column index
    col_positions = [header_vals.index(col) for col in needed]

    # 5) Fix merged‐cell misalignment:
    #    if the “header” column j is blank for all data rows
    #    but j+1 has real values, shift to j+1
    for idx, j in enumerate(col_positions):
        if needed[idx] == "Bacteria Species":
            # An intentionally blank organism column is not a shifted barcode
            # or another metadata column. Keep its missing-input provenance.
            continue
        block = raw.iloc[data_start : data_start + 10, j]
        # is it entirely blank?
        if all(pd.isna(x) or str(x).strip() == "" for x in block):
            # does the next column have any data?
            if j + 1 < n_cols and any(
                not (pd.isna(x) or str(x).strip() == "")
                for x in raw.iloc[data_start : data_start + 10, j + 1]
            ):
                col_positions[idx] = j + 1

    # 6) Slice out exactly those columns, build a DataFrame
    data_block = raw.iloc[data_start:, col_positions]
    df = pd.DataFrame(data_block.values, columns=needed)

    # 7) Clean up every cell
    df = df.apply(lambda col: col.map(normalize))

    # The controlled workbook contains reagent/protocol sections below the
    # sample table. A blank sample row ends that table; dropping blanks first
    # would join those unrelated sections back into the sample list. Barcode
    # presence must not determine this boundary (uploaded samples may omit it).
    table_started = False
    for position, (_, row) in enumerate(df.iterrows()):
        has_sample_fields = bool(row["Accession Number"] or row["Bacteria Species"])
        if table_started and not has_sample_fields:
            df = df.iloc[:position].copy()
            break
        if has_sample_fields:
            table_started = True

    if "Other Bacteria Species" in df.columns:
        other_mask = df["Bacteria Species"].str.casefold() == "other"
        supplied_other = other_mask & df["Other Bacteria Species"].astype(bool)
        df.loc[supplied_other, "Bacteria Species"] = df.loc[
            supplied_other, "Other Bacteria Species"
        ]

    # 8) Illumina rows require a valid well. Nanopore barcode is optional for
    # explicitly sample-mapped uploads; when absent the sample folder is used.
    if plat == "illumina":
        df = df[df[idx_col].astype(bool) & df[idx_col].str.match(pattern)]

    # 9) Drop any row missing the two critical fields
    df = df[df["Accession Number"].astype(bool)]

    # 10) Return a list of TAB‐joined strings

    if plat == "nanopore":
        out = df.apply(
            lambda row: (
                f"{row['Accession Number']}\t{row['Bacteria Species']}\t"
                f"{row[idx_col] if idx_col in df.columns and row[idx_col] else row['Accession Number']}"
            ),
            axis=1
        ).tolist()
    elif plat == "illumina":
        out = df.apply(
            lambda row: f"{row['Accession Number']}\t{row['Bacteria Species']}",
            axis=1
        ).tolist()
    
    return out

def download_s3(bucket: str, prefix: str, target_root: str) -> str:
    s3 = boto3.client("s3")
    norm_prefix = prefix.strip("/")
    base_name = norm_prefix.split("/")[-1]
    local_base = os.path.join(target_root, base_name)
    os.makedirs(local_base, exist_ok=True)

    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=norm_prefix + "/"):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.endswith("/"):
                continue
            rel_path = key[len(norm_prefix):].lstrip("/")
            if not rel_path:
                continue
            dest_path = os.path.join(local_base, rel_path)
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)
            s3.download_file(bucket, key, dest_path)

    return local_base


init(autoreset=True)

class TqdmMinutes(tqdm):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._rate_format = "{rate:.2f} min/genome"

    def format_dict(self):
        d = super().format_dict
        rate_min = (self.format_dict['elapsed'] / self.n) / 60 if self.n else 0
        d['rate'] = rate_min
        return d

    def format_meter(self, n, total, elapsed, ncols=None):
        rate_min = (elapsed / n) / 60 if n else 0
        return super().format_meter(n, total, elapsed, ncols).replace("it/s", f"{rate_min:.2f} min/genome")
    
# Compress files at the end of run
def compress_qc_fastqs(
    qc_out: str,
    genomes_dir: str,
    cpus_per_sample: int = 4,
    require_assembly: bool = True,
    recursive: bool = False,
    remove_source: bool = True,
    logger_fn: Optional[Callable[[str, str], None]] = None,
) -> Dict[str, int]:
    """
    Compress uncompressed FASTQ files under qc_out/<sample> to .gz in parallel (threads).

    - Uses pigz if available (multi-threaded), else gzip.
    - Atomic write: *.gz.tmp -> os.replace(...).
    - Verifies integrity (pigz/gzip -t) before deleting source.
    - Avoids CPU oversubscription: workers ~= total_cpus // threads_per_file.
    - Logs via standard `logging` so records go into your main log file.
      Optionally accepts `logger_fn(message, level)`; if provided, it’s used first.

    Returns: {"total": N, "ok": K, "failed": F, "skipped": S}
    """
    log = logging.getLogger(__name__)

    def _log(msg: str, level: str = "Norm") -> None:
        if logger_fn:
            try:
                logger_fn(msg, level)   # e.g., compat_logger(..., message_type=level)
                return
            except TypeError:
                try:
                    logger_fn(msg)      # accept simple print-like callables too
                    return
                except Exception:
                    pass
        log.log(_LEVEL_MAP.get(level, logging.INFO), msg)

    compressor = shutil.which("pigz") or shutil.which("gzip")
    if not compressor:
        _log("Neither pigz nor gzip found on PATH; skipping compression.", "Warn")
        return {"total": 0, "ok": 0, "failed": 0, "skipped": 0}

    # ---- Gather jobs ----------------------
    jobs: List[Tuple[str, str]] = []  # (in_fastq, out_gz)
    skipped = 0

    for sample in sorted(os.listdir(qc_out)):
        sample_dir = os.path.join(qc_out, sample)
        if not os.path.isdir(sample_dir):
            continue

        if require_assembly:
            asm_path = os.path.join(genomes_dir, f"{sample}.fasta")
            if not (os.path.exists(asm_path) and os.path.getsize(asm_path) > 0):
                skipped += 1
                continue

        if recursive:
            for root, _, files in os.walk(sample_dir):
                for name in files:
                    if name.endswith(".fastq"):         # add ".fq" here if you want
                        in_path = os.path.join(root, name)
                        out_gz = in_path + ".gz"
                        if not os.path.exists(out_gz):
                            jobs.append((in_path, out_gz))
        else:
            for name in os.listdir(sample_dir):
                if name.endswith(".fastq"):
                    in_path = os.path.join(sample_dir, name)
                    out_gz = in_path + ".gz"
                    if not os.path.exists(out_gz):
                        jobs.append((in_path, out_gz))

    if not jobs:
        _log("No uncompressed FASTQ files found to compress.", "Norm")
        return {"total": 0, "ok": 0, "failed": 0, "skipped": skipped}

    # ---- Parallelism budgeting --------------------------------------------
    try:
        threads_per_file = max(1, int(cpus_per_sample))
    except Exception:
        threads_per_file = 4
    if os.path.basename(compressor) != "pigz":
        threads_per_file = 1  # gzip is single-threaded

    total_cpus = os.cpu_count() or 1
    max_workers = max(1, min(len(jobs), max(1, total_cpus // threads_per_file)))

    _log(
        f"Compressing {len(jobs)} FASTQ files with {os.path.basename(compressor)} "
        f"(workers={max_workers}, threads/file={threads_per_file}).",
        "Norm",
    )

    # ---- Helpers -----------------------------------------------------------
    def _integrity_ok(gz_path: str) -> bool:
        tester = shutil.which("pigz") or shutil.which("gzip")
        if not tester:
            return True
        r = subprocess.run(
            [tester, "-t", gz_path],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return r.returncode == 0

    def _compress_one(in_path: str, out_gz: str) -> str:
        gz_tmp = out_gz + ".tmp"
        # Remove stale tmp if present
        try:
            if os.path.exists(gz_tmp):
                os.remove(gz_tmp)
        except Exception:
            pass

        cmd = [compressor]
        if os.path.basename(compressor) == "pigz" and threads_per_file > 1:
            cmd += ["-p", str(threads_per_file)]
        cmd += ["-c", in_path]

        # Primary attempt
        try:
            with open(gz_tmp, "wb") as gz_out:
                subprocess.run(
                    cmd,
                    check=True,
                    stdin=subprocess.DEVNULL,
                    stdout=gz_out,
                    stderr=subprocess.DEVNULL,
                )
        except subprocess.CalledProcessError:
            # Fallback to gzip if pigz failed
            if os.path.basename(compressor) == "pigz":
                gz = shutil.which("gzip")
                if not gz:
                    raise
                with open(gz_tmp, "wb") as gz_out:
                    subprocess.run(
                        [gz, "-c", in_path],
                        check=True,
                        stdin=subprocess.DEVNULL,
                        stdout=gz_out,
                        stderr=subprocess.DEVNULL,
                    )
            else:
                raise

        # Verify gzip integrity
        if not _integrity_ok(gz_tmp):
            try:
                os.remove(gz_tmp)
            except FileNotFoundError:
                pass
            raise RuntimeError(f"Compression test failed for {in_path}")

        os.replace(gz_tmp, out_gz)  # atomic publish

        if remove_source:
            try:
                os.remove(in_path)
            except FileNotFoundError:
                pass

        return out_gz

    # ---- Run in thread pool -----------------------------------------------
    ok = 0
    failed = 0
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = {ex.submit(_compress_one, src, dst): (src, dst) for (src, dst) in jobs}
        for i, fut in enumerate(as_completed(futs), start=1):
            src, _dst = futs[fut]
            try:
                fut.result()
                ok += 1
                if i % 10 == 0 or i == len(jobs):
                    _log(f"Compressed {ok}/{len(jobs)} files.", "Norm")
            except Exception as e:
                failed += 1
                _log(f"Compression failed for {src}: {e}", "Warn")

    if failed:
        _log(f"End-of-pipeline compression finished with {failed} failure(s).", "Warn")
    else:
        _log("End-of-pipeline compression finished successfully.", "Pass")

    return {"total": len(jobs), "ok": ok, "failed": failed, "skipped": skipped}

# ---------- fast, safe open for SeqIO (supports .gz via pigz/gzip -dc) ----------
def _open_text_fastq_for_seqio(path: str, pigz_threads: int = 1) -> Tuple[io.TextIOBase, Optional[subprocess.Popen]]:
    """
    Returns (text_handle, proc). If proc is not None, caller should wait() it after use.
    """
    proc = None
    if path.endswith(".gz"):
        dec = shutil.which("pigz") or shutil.which("gzip")
        if dec:
            # robust int for pigz threads
            try:
                t = max(1, int(pigz_threads))
            except Exception:
                t = 1
            cmd = [dec]
            if os.path.basename(dec) == "pigz" and t > 1:
                cmd += ["-p", str(t)]
            cmd += ["-dc", path]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            assert proc.stdout is not None
            fh = io.TextIOWrapper(io.BufferedReader(proc.stdout, 1024 * 1024),
                                  encoding="ascii", newline="")
            return fh, proc
        else:
            fh = io.TextIOWrapper(gzip.open(path, "rb"), encoding="ascii", newline="")
            return fh, None
    else:
        return open(path, "rt", encoding="ascii", newline=""), None


# =========================== API 1: fastq_metrics ============================

def fastq_metrics(path: str, use_external_decompress: bool = True, pigz_threads=1):
    """
    Streaming metrics from FASTQ/FASTQ.GZ using Biopython/SeqIO.
    Returns a dict with:
      reads, total_bases,
      mean_q_read_np   -> NanoPlot/MinKNOW-style mean read quality (per-read prob avg -> Phred)
      mean_q_base_np   -> Phred of base-weighted mean error prob across all bases
      mean_q_per_base_arith -> arithmetic mean of Phred across bases (legacy/comparison)
      mean_q_per_read_arith -> arithmetic mean of per-read mean Phred (legacy/comparison)
    """
    if path.endswith(".gz") and use_external_decompress:
        fh, proc = _open_text_fastq_for_seqio(path, pigz_threads=pigz_threads)
    else:
        fh, proc = open(path, "rt", encoding="ascii", newline=""), None

    reads = 0
    total_bases = 0

    # NanoPlot/MinKNOW accumulators
    p_read_sum = 0.0           # sum over reads of mean per-base error prob
    p_base_sum = 0.0           # sum over all bases of error prob (i.e., Σ_read p_read * L)

    # Arithmetic means (for comparison with FastQC-like summaries)
    qsum_all = 0               # Σ over all bases of Phred
    read_mean_qsum = 0.0       # Σ over reads of (mean Phred per read)

    try:
        for rec in SeqIO.parse(fh, "fastq"):   # PHRED+33 integers
            quals = rec.letter_annotations.get("phred_quality")
            if not quals:
                continue
            L = len(quals)
            if L <= 0:
                continue

            # NanoPlot/MinKNOW-style per-read mean error prob
            p_read = sum(10 ** (-q / 10.0) for q in quals) / L
            p_read_sum += p_read
            p_base_sum += p_read * L

            # Arithmetic accumulators
            qsum = sum(quals)
            qsum_all += qsum
            read_mean_qsum += (qsum / L)

            reads += 1
            total_bases += L
    finally:
        fh.close()
        if proc is not None:
            try: proc.wait()
            except Exception: pass

    # Convert prob means back to Phred (guard against zeros)
    mean_q_read_np = (-10.0 * math.log10(p_read_sum / reads)) if reads and p_read_sum > 0 else 0.0
    mean_q_base_np = (-10.0 * math.log10(p_base_sum / total_bases)) if total_bases and p_base_sum > 0 else 0.0

    # Arithmetic means
    mean_q_per_base_arith = (qsum_all / total_bases) if total_bases else 0.0
    mean_q_per_read_arith = (read_mean_qsum / reads) if reads else 0.0

    return {
        "reads": reads,
        "total_bases": total_bases,
        "mean_q_read_np": mean_q_read_np,
        "mean_q_base_np": mean_q_base_np,
        "mean_q_per_base_arith": mean_q_per_base_arith,
        "mean_q_per_read_arith": mean_q_per_read_arith,
    }


# ================= API 2: run_filtlong_and_write_with_metrics ================

def filtlong_with_metrics(filtlong_cmd, out_fastq_path: str):
    """
    Run filtlong, stream its FASTQ output, write UNCOMPRESSED out_fastq_path,
    and compute NanoPlot/MinKNOW-style metrics in one pass.
    Returns the same keys as fastq_metrics(...).
    """
    p = subprocess.Popen(
        filtlong_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,       # drain stderr so it can't block
        stdin=subprocess.DEVNULL,
        bufsize=1024 * 1024,
        text=False,
    )
    assert p.stdout is not None

    # drain stderr in background
    def _drain_err(stream):
        while True:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break

    t = threading.Thread(target=_drain_err, args=(p.stderr,), daemon=True)
    t.start()

    # Wrap stdout for text parsing
    fh = io.TextIOWrapper(io.BufferedReader(p.stdout, 1024 * 1024),
                          encoding="ascii", newline="")

    reads = 0
    total_bases = 0
    p_read_sum = 0.0
    p_base_sum = 0.0
    qsum_all = 0
    read_mean_qsum = 0.0

    with open(out_fastq_path, "wt", encoding="ascii", newline="") as out:
        for rec in SeqIO.parse(fh, "fastq"):
            quals = rec.letter_annotations.get("phred_quality")
            if not quals:
                continue
            L = len(quals)
            if L <= 0:
                continue

            # NanoPlot-style accumulators
            p_read = sum(10 ** (-q / 10.0) for q in quals) / L
            p_read_sum += p_read
            p_base_sum += p_read * L

            # Arithmetic accumulators (optional, for comparison)
            qsum = sum(quals)
            qsum_all += qsum
            read_mean_qsum += (qsum / L)

            reads += 1
            total_bases += L

            # Write normalized 4-line FASTQ
            SeqIO.write(rec, out, "fastq")

    ret = p.wait()
    t.join()
    if ret != 0:
        raise subprocess.CalledProcessError(ret, filtlong_cmd)

    mean_q_read_np = (-10.0 * math.log10(p_read_sum / reads)) if reads and p_read_sum > 0 else 0.0
    mean_q_base_np = (-10.0 * math.log10(p_base_sum / total_bases)) if total_bases and p_base_sum > 0 else 0.0
    mean_q_per_base_arith = (qsum_all / total_bases) if total_bases else 0.0
    mean_q_per_read_arith = (read_mean_qsum / reads) if reads else 0.0

    return {
        "reads": reads,
        "total_bases": total_bases,
        "mean_q_read_np": mean_q_read_np,
        "mean_q_base_np": mean_q_base_np,
        "mean_q_per_base_arith": mean_q_per_base_arith,
        "mean_q_per_read_arith": mean_q_per_read_arith,
    }


## For merging files
def merge_gz_fastqs(raw_folder: str, output_gz: str, *, verify: bool = False) -> str:
    """
    Merge all *.fastq.gz/.fq.gz files in deterministic natural-sort order.
    Writes to <output_gz>.tmp then atomically renames. Returns output path.
    """
    if not output_gz.endswith(".gz"):
        raise ValueError("output_gz must end with .gz")

    out_dir = os.path.dirname(os.path.abspath(output_gz)) or "."
    os.makedirs(out_dir, exist_ok=True)

    tmp_out = output_gz + ".tmp"
    inputs = [
        os.path.join(raw_folder, name)
        for name in os.listdir(raw_folder)
        if name.endswith((".fastq.gz", ".fq.gz"))
    ]
    inputs.sort(key=_natural_key)
    if not inputs:
        raise FileNotFoundError(f"No .fastq.gz/.fq.gz files found in {raw_folder}")

    with open(tmp_out, "wb") as out_fh:
        for in_path in inputs:
            with open(in_path, "rb") as in_fh:
                shutil.copyfileobj(in_fh, out_fh, length=8 * 1024 * 1024)

    if verify:
        gz = shutil.which("pigz") or shutil.which("gzip")
        if gz:
            subprocess.run([gz, "-t", tmp_out], check=True)

    os.replace(tmp_out, output_gz)
    return output_gz

#S3 sources: merge_gz_fastqs_s3

def _natural_key(value: str):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", value)]


def _gz_integrity_ok(path: str) -> bool:
    gz = shutil.which("pigz") or shutil.which("gzip")
    if not gz:
        return True
    return subprocess.run(
        [gz, "-t", path],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    ).returncode == 0

def _parse_s3_uri(uri: str) -> Tuple[str, str]:
    if not uri.startswith("s3://"):
        raise ValueError(f"Not an s3 URI: {uri}")
    without = uri[5:]
    bucket, _, key = without.partition("/")
    if not bucket:
        raise ValueError(f"Invalid s3 URI: {uri}")
    return bucket, key

def _list_s3_fastq_gz(session, bucket: str, prefix: str) -> List[Tuple[str, str]]:
    s3 = session.client("s3")
    paginator = s3.get_paginator("list_objects_v2")
    keys: List[Tuple[str, str]] = []
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            k = obj["Key"]
            if k.endswith((".fastq.gz", ".fq.gz")):
                keys.append((bucket, k))
    keys.sort(key=lambda bk: _natural_key(bk[1]))
    return keys

def merge_gz_fastqs_s3(
    s3_sources_or_prefix: Iterable[str] | str,
    output_path: str,
    aws_profile: Optional[str] = None,
    region: Optional[str] = None,
    verify: bool = False,
    logger_fn: Optional[Callable[[str, str], None]] = None,
) -> dict:
    """
    Fast merge of S3 .fastq.gz/.fq.gz objects -> local .gz or uncompressed FASTQ.
    Requires boto3.
    """
    def log(msg: str, level: str = "Norm"):
        if logger_fn:
            try:
                logger_fn(msg, level)
            except Exception:
                pass

    try:
        import boto3
    except ImportError:
        raise RuntimeError("boto3 is required for merge_gz_fastqs_s3.")

    session = boto3.Session(profile_name=aws_profile, region_name=region)

    # Build (bucket, key) list
    inputs: List[Tuple[str, str]] = []
    if isinstance(s3_sources_or_prefix, str):
        bucket, key = _parse_s3_uri(s3_sources_or_prefix)
        if key.endswith("/"):
            inputs = _list_s3_fastq_gz(session, bucket, key)
        else:
            inputs = [(bucket, key)]
    else:
        for uri in s3_sources_or_prefix:
            b, k = _parse_s3_uri(str(uri))
            inputs.append((b, k))
        inputs.sort(key=lambda bk: _natural_key(bk[1]))

    if not inputs:
        raise FileNotFoundError("No S3 .fastq.gz/.fq.gz objects found to merge.")

    s3 = session.client("s3")
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    tmp_out = output_path + ".tmp"

    try:
        if output_path.endswith(".gz"):
            log(f"Merging {len(inputs)} S3 gz FASTQs -> {output_path} (raw concat).")
            total = 0
            with open(tmp_out, "wb") as out:
                for bucket, key in inputs:
                    obj = s3.get_object(Bucket=bucket, Key=key)
                    body = obj["Body"]
                    for chunk in iter(lambda: body.read(8 * 1024 * 1024), b""):
                        out.write(chunk)
                        total += len(chunk)
            if verify and not _gz_integrity_ok(tmp_out):
                raise RuntimeError("Gzip integrity test failed for merged output.")
            os.replace(tmp_out, output_path)
            return {"inputs": len(inputs), "bytes_written": total, "output": output_path}
        else:
            log(f"Merging {len(inputs)} S3 gz FASTQs -> {output_path} (decompress).")
            total = 0
            import gzip  # python gzip to decompress stream
            with open(tmp_out, "wb") as out:
                for bucket, key in inputs:
                    obj = s3.get_object(Bucket=bucket, Key=key)
                    body = obj["Body"]
                    # Wrap the streaming body in a file-like for gzip
                    with gzip.GzipFile(fileobj=body, mode="rb") as gz:
                        shutil.copyfileobj(gz, out, length=8 * 1024 * 1024)
            total = os.path.getsize(tmp_out)
            os.replace(tmp_out, output_path)
            return {"inputs": len(inputs), "bytes_written": total, "output": output_path}
    finally:
        try:
            if os.path.exists(tmp_out):
                os.remove(tmp_out)
        except Exception:
            pass

##########TYPING ANALYSIS HELPERS ##########

def _fmt_int_commas(val: str) -> str:
    """Format integers with thousands separators; pass through NA/blank."""
    v = (val or "").strip()
    if v == "" or v.upper() == "NA":
        return v or "NA"
    try:
        # SKA 'Distance' is integer; guard against accidental floats/strings
        i = int(float(v))
        return f"{i:,}"
    except Exception:
        return v

# Typer reports
# bactipipe/scripts/utils.py  (append these helpers)
# If your runner is available, use it; otherwise fallback to subprocess
try:
    from .runner import run as _run  # signature: (cmd: List[str], cwd: Optional[Path], log, env_name: Optional[str]) -> (rc, out)
except Exception:
    _run = None  # fallback below


def _run_cmd(cmd: List[str], env_name: Optional[str], logger) -> Tuple[int, str]:
    """
    Run a command, optionally inside a conda/mamba env via 'conda run -n <env>'.
    Uses .runner.run if available; else a minimal subprocess fallback.
    """
    if _run is not None:
        return _run(cmd, cwd=None, log=logger, env_name=env_name)

    # minimal fallback (mirrors your _compose_tool_cmd)
    def _compose(cmd_: List[str], env: Optional[str]) -> List[str]:
        if not env:
            return cmd_
        for tool in ("conda", "micromamba", "mamba"):
            if shutil.which(tool):
                if tool == "conda":
                    return [tool, "run", "--no-capture-output", "-n", env] + cmd_
                return [tool, "run", "-n", env] + cmd_
        return cmd_
    import shutil
    full = _compose(cmd, env_name)
    try:
        p = subprocess.run(full, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        return p.returncode, p.stdout
    except Exception as e:
        return 127, f"[runner] Failed to execute: {e}"


# ---------------------------
# Version normalization
# ---------------------------

_VERSION_RE = re.compile(r"(?:^|[^A-Za-z0-9_])(?P<v>v?\d+(?:\.\d+){0,3})(?:[^A-Za-z0-9_.]|$)")

def normalize_version_string(raw: str) -> str:
    """
    Extracts the first version-ish token and returns it prefixed with 'v'.
    Examples:
      'SeqSero2_package.py 1.2.1'   -> 'v1.2.1'
      'Kleborate v3.2.4'            -> 'v3.2.4'
      'skani 0.3.0'                 -> 'v0.3.0'
      'v1.2.0'                      -> 'v1.2.0'
    """
    if not raw:
        return "v?"
    m = _VERSION_RE.search(raw.strip())
    if not m:
        return "v?"
    token = m.group("v")
    return token if token.startswith("v") else f"v{token}"


def get_tool_version(
    exe: str,
    logger,
    env_name: Optional[str] = None,
    flags: Optional[List[str]] = None,
    alt_flags: Optional[List[str]] = None,
    manual: Optional[str] = None,
) -> str:
    """
    Tries flags then alt_flags; on success normalizes to 'vX.Y.Z'.
    On failure returns 'manual' if provided, else 'v?'.
    """
    flags = flags or ["--version"]
    alt_flags = alt_flags or []
    tried: List[List[str]] = []
    for fl in (flags, alt_flags):
        if not fl:
            continue
        tried.append(fl)
        rc, out = _run_cmd([exe] + fl, env_name, logger)
        if rc == 0 and out and out.strip():
            return normalize_version_string(out.splitlines()[0])
    if manual:
        return manual if manual.startswith("v") else f"v{manual}"
    logger.debug(f"[version] Failed to detect version for {exe} with flags {tried}; using v?")
    return "v?"

# ---------------------------
# Tools used & version collection
# ---------------------------

def get_python_version() -> str:
    import sys
    v = sys.version_info
    return f"v{v.major}.{v.minor}.{v.micro}"

def get_biopython_version(logger=None) -> str:
    try:
        import Bio  # type: ignore
        ver = getattr(Bio, "__version__", None) or "N/A"
        return ver if ver.startswith("v") else f"v{ver}"
    except Exception as e:
        if logger:
            logger.debug(f"Biopython not importable: {e}")
        return "N/A"

def detect_used_tools(
    paths: Dict[str, Path],
    logger,
    *,
    env_overrides: Optional[Dict[str, Optional[str]]] = None,
) -> Dict[str, Dict[str, Optional[str]]]:
    """
    Build a 'used' mapping: {ToolName -> {exe, env} or {version}} based on outputs in `paths`.
    - Always includes Python and Biopython (prefilled versions).
    - Includes MLST and SerotypeFinder as manual versions if used.
    - Detects skani for both ref-based ANI and triangle runs.
    """
    env_overrides = env_overrides or {}
    used: Dict[str, Dict[str, Optional[str]]] = {}

    # Always add Python & Biopython (prefilled versions)
    # used["Python"] = {"exe": None, "env": None, "version": get_python_version()}
    used["Biopython"] = {"exe": None, "env": None, "version": get_biopython_version(logger)}

    # SeqSero2 present if serotype/…/seqsero2 outputs exist
    s_root = paths.get("serotype_root")
    if s_root and s_root.exists():
        any_seqsero = any(
            (p / "seqsero2" / "SeqSero_result.tsv").exists()
            or (p / "seqsero2" / "results_tab.tsv").exists()
            or (p / "seqsero2" / "seqsero2.stdout.txt").exists()
            for p in s_root.iterdir() if p.is_dir()
        )
        if any_seqsero:
            used["SeqSero2"] = {"exe": "SeqSero2_package.py", "env": env_overrides.get("seqsero2")}

    # cgMLSTFinder if cgmlst TSV or per-sample summaries exist
    cgmlst_tsv = paths.get("cgmlst_tsv")
    cgmlst_root = paths.get("cgmlst_root")
    if (cgmlst_tsv and cgmlst_tsv.exists()) or (
        cgmlst_root and cgmlst_root.exists() and any(
            (p / "salmonella_summary.txt").exists() or (p / "cgmlst.stdout.txt").exists()
            for p in cgmlst_root.iterdir() if p.is_dir()
        )
    ):
        used["cgMLSTFinder"] = {"exe": "cgMLST.py", "env": env_overrides.get("cge")}

    # MLST if mlst TSV or per-sample results.txt exist (manual version pin)
    mlst_tsv = paths.get("mlst_tsv")
    mlst_root = paths.get("mlst_root")
    if (mlst_tsv and mlst_tsv.exists()) or (
        mlst_root and mlst_root.exists() and any(
            (p / "results.txt").exists() for p in mlst_root.iterdir() if p.is_dir()
        )
    ):
        used["MLST"] = {"exe": None, "env": None, "version": "v2.0.9"}

    # SerotypeFinder if any serotypefinder results exist (manual version pin)
    if s_root and s_root.exists():
        if any((p / "serotypefinder" / "results_tab.tsv").exists()
               for p in s_root.iterdir() if p.is_dir()):
            used["SerotypeFinder"] = {"exe": None, "env": None, "version": "v2.0.1"}

    # Kleborate if any kleborate dir exists
    if s_root and s_root.exists():
        if any((p / "kleborate").exists() for p in s_root.iterdir() if p.is_dir()):
            used["Kleborate"] = {"exe": "kleborate", "env": env_overrides.get("kleborate")}

    # SKA2 if distances TSV exists
    ska_dist = paths.get("ska_distances")
    if ska_dist and ska_dist.exists():
        used["SKA2"] = {"exe": "ska", "env": env_overrides.get("ska")}

    # skani: present if ref ANI table exists OR triangle pairs exist
    ani_present = bool(paths.get("ani_tsv") and paths["ani_tsv"].exists())
    tri_present = any(paths.get(k) and paths[k].exists()
                      for k in ("skani_pairs", "skani_pairs_samples"))
    if ani_present or tri_present:
        used["skani"] = {"exe": "skani", "env": env_overrides.get("skani")}

    return used

def _norm_ver(v: Optional[str]) -> str:
    v = (v or "").strip()
    if not v:
        return "N/A"
    return v if v.startswith("v") else f"v{v}"

def collect_tool_versions(used: Dict[str, Dict[str, Optional[str]]], logger) -> Dict[str, str]:
    """
    Convert a 'used' mapping -> {name: version}.
    Supports entries with a prefilled 'version' OR with ('exe','env') to probe via get_tool_version.
    Normalizes everything to start with 'v', returns 'N/A' on failure.
    """
    versions: Dict[str, str] = {
        "BactiPipe": bactipipe_version,
        "BactiPipe database": database_identity_label(),
    }

    for name, meta in used.items():
        # Be tolerant: someone may accidentally pass a string
        if not isinstance(meta, dict):
            versions[name] = _norm_ver(str(meta))
            continue

        # Prefilled version takes precedence (Python, Biopython, MLST, SerotypeFinder, etc.)
        pref = meta.get("version")
        if pref:
            versions[name] = _norm_ver(pref)
            continue

        exe = meta.get("exe")
        env = meta.get("env")

        # Manual pins (if you prefer to keep them here instead of prefill in 'used')
        if name == "MLST":
            versions[name] = "v2.0.9"
            continue
        if name == "SerotypeFinder":
            versions[name] = "v2.0.1"
            continue

        if not exe:
            versions[name] = "N/A"
            continue

        # Probe
        try:
            # Your utils.get_tool_version already normalizes, but we normalize again defensively
            v = get_tool_version(exe, logger, env_name=env)  # <- keep existing signature
            versions[name] = _norm_ver(v)
        except Exception:
            versions[name] = "N/A"

    return versions

# ---------------------------
# Report header + PDF
# ---------------------------

def make_report_header_block(
    *,
    title: str,
    run_name: str,
    tech_name: str,
    accession: str = "",
) -> str:
    """
    Single-block header (accession inline). Matches your request.
    """
    date_str = datetime.now().strftime("%Y-%m-%d %H:%M")
    hostname = socket.gethostname()
    try:
        ip_addr = socket.gethostbyname(hostname)
    except Exception:
        ip_addr = "unknown"

    header = f"""
=======================================
Penn State Animal Diagnostic Laboratory
BactiPipe - {title}
Run name: {run_name}
Accession: {accession or '-'}
Date: {date_str}
Analysis by: {tech_name}
Server Host Name: {hostname}
=======================================
""".strip("\n")
    return header

def render_pdf_type_genomes_full(
    *,
    final_tsv: Path,
    out_pdf: Path,
    title: str,
    header_text: str,
    tool_versions: Dict[str, str],
    run_name: str,
    accession: str,
) -> None:
    """
    Minimal, robust renderer:
      - No reordering of columns (uses whatever is in the TSV)
      - Wraps only 'Specimen' column
      - Centers numeric columns: Size (Mb), GC (%), MLST, cgMLST, SNP distance, ANI (%)
      - Forces 'SNP distance' header to two lines (SNP<br/>distance)
      - Fits every table within page width (landscape Letter)
    """
    if not final_tsv.exists():
        raise FileNotFoundError(f"Final TSV not found: {final_tsv}")
    # ---------- Parse final.tsv into sections (## Title lines start sections) ----------
    raw_lines = final_tsv.read_text().splitlines()
    sections: List[Tuple[str, List[List[str]]]] = []
    current_title = "Strain Relatedness Summary"
    current_lines: List[str] = []

    def _flush_current():
        nonlocal current_title, current_lines
        if current_lines:
            rows: List[List[str]] = []
            rdr = csv.reader(current_lines, delimiter="\t")
            for r in rdr:
                if r and any((c or "").strip() for c in r):
                    rows.append([c.strip() for c in r])
            if rows:
                sections.append((current_title, rows))
        current_lines = []

    for line in raw_lines:
        if line.startswith("## "):
            _flush_current()
            current_title = line[3:].strip() or "Section"
        else:
            if line.strip() == "" and not current_lines:
                continue
            current_lines.append(line)
    _flush_current()
    if not sections:
        raise RuntimeError("No sectioned tables detected in final TSV.")

    # ---------- Styles ----------
    styles = getSampleStyleSheet()

    title_style = ParagraphStyle(
        "TitleCentered",
        parent=styles["Heading1"],
        alignment=TA_CENTER,
        fontName="Helvetica-Bold",
        fontSize=16,
        spaceAfter=6,
    )
    hdr_title = ParagraphStyle(
        "HeaderTitle",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=10,
        leading=12,
    )
    hdr_block = ParagraphStyle(
        "HeaderBlock",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=9,
        leading=11,
    )
    hdr_sep = ParagraphStyle(
        "HeaderSep",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=9,
        leading=11,
    )
    section_h = ParagraphStyle(
        "SectionHeading",
        parent=styles["Heading2"],
        fontName="Helvetica-Bold",
        fontSize=12,
        spaceBefore=8,
        spaceAfter=4,
    )
    header_cell = ParagraphStyle(
        "HeaderCell",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=9,
        leading=11,
        alignment=TA_CENTER,  # header text centered
    )
    cell_p = ParagraphStyle(
        "Cell",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=9,
        leading=11,
        wordWrap="LTR",
        alignment=TA_LEFT,  # body defaults to left
    )
    cell_center = ParagraphStyle(
        "CellCenter",
        parent=cell_p,
        alignment=TA_CENTER,
    )
    # Header cell styles
    header_cell_left = ParagraphStyle(
        "HeaderCellLeft",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=9,
        leading=11,
        alignment=TA_LEFT,
    )
    header_cell_center = ParagraphStyle(
        "HeaderCellCenter",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=9,
        leading=11,
        alignment=TA_CENTER,
    )

    # ---------- Footer (right-aligned) ----------
    def _on_page(canvas, doc):
        canvas.setFont("Helvetica", 9)
        footer = f"BactiPipe - Relatedness | {run_name}/{accession}  | Page {doc.page}"
        x = doc.pagesize[0] - doc.rightMargin
        canvas.drawRightString(x, 0.45 * inch, footer)

    # ---------- Header band helpers ----------
    def _fit_sep(panel_width: float, font="Helvetica", size=9, scale=0.50,
                 min_chars=4, max_chars=60) -> str:
        usable = max(1.0, panel_width - 10)
        per_char = max(1.0, stringWidth("=", font, size))
        n = int((usable / per_char) * scale)
        n = max(min_chars, min(max_chars, n))
        return "=" * n

    # ---------- Column width calculator (kept compact) ----------
    def _cell_plain_text(x) -> str:
        if isinstance(x, Paragraph):
            return x.getPlainText()
        return "" if x is None else str(x)

    def _column_widths(data: List[List[object]], max_width: float) -> List[float]:
        if not data or not data[0]:
            return []
        ncols = len(data[0])
        widths = [0.0] * ncols
        # Estimate widths from up to first 200 body rows + header
        sample_rows = data[:min(len(data), 200)]
        for row in ([data[0]] + sample_rows):
            for j, cell in enumerate(row):
                s = _cell_plain_text(cell)
                w = stringWidth(s, "Helvetica", 9) + 12  # padding
                widths[j] = max(widths[j], w)

        # Enforce semantic minimums
        header_txts = [re.sub(r"<.*?>", "", _cell_plain_text(h)).strip() for h in data[0]]
        lower = [h.lower() for h in header_txts]

        def ensure_min(col_name: str, min_inch: float):
            try:
                i = lower.index(col_name.lower())
                widths[i] = max(widths[i], min_inch * inch)
            except ValueError:
                pass

        # Make Specimen generous; numeric cols moderate
        ensure_min("Sample ID", 1.30)
        ensure_min("Isolate ID", 1.30)
        ensure_min("Specimen", 2.0)
        ensure_min("Serotype", 1.25)
        ensure_min("SNP distance", 1.20)
        ensure_min("ANI (%)", 1.00)
        ensure_min("Size (Mb)", 1.00)
        ensure_min("GC (%)", 1.00)
        ensure_min("MLST", 0.8)
        ensure_min("cgMLST", 0.9)

        total = sum(widths) or 1.0
        if total <= max_width:
            return widths

        # Scale down proportionally to fit within max_width
        scale = max_width / total
        widths = [w * scale for w in widths]

        # If tiny overflow due to rounding, shave proportionally (respect a minimal floor)
        floor = 0.60 * inch
        total = sum(widths)
        if total > max_width:
            over = total - max_width
            flex = [max(0, w - floor) for w in widths]
            flex_total = sum(flex) or 1.0
            widths = [max(floor, w - over * (max(0, w - floor) / flex_total)) for w in widths]

        # As an extra guard: cap any single column to 45% of the table width
        cap = max_width * 0.45
        clipped = False
        for i, w in enumerate(widths):
            if w > cap:
                widths[i] = cap
                clipped = True
        if clipped:
            # re-distribute remainder lightly (optional; can be omitted to keep it simple)
            pass

        return widths

    # ---------- Build the document ----------
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(
        str(out_pdf),
        pagesize=landscape(letter),
        leftMargin=0.6 * inch,
        rightMargin=0.6 * inch,
        topMargin=0.6 * inch,
        bottomMargin=0.6 * inch,
    )
    story: List = []

    # Title
    story.append(Paragraph(title, title_style))
    story.append(Spacer(1, 6))

    # Two-block header (left: Run Metadata; right: Tool versions)
    col_widths_hdr = [doc.width * 0.56, doc.width * 0.44]

    # LEFT
    left_flow: List[Paragraph] = [Paragraph("Run Metadata", hdr_title)]
    for ln in (header_text or "").strip("\n").splitlines():
        if re.fullmatch(r"=+", ln.strip()):
            left_flow.append(Paragraph(ln.strip(), hdr_sep))
        else:
            left_flow.append(Paragraph(str(ln), hdr_block))

    # RIGHT
    right_sep = _fit_sep(col_widths_hdr[1])
    right_flow: List[Paragraph] = [Paragraph("Software Version Numbers", hdr_title),
                                   Paragraph(right_sep, hdr_sep)]
    # Use order given by tool_versions (dict insertion order)
    for name, ver in tool_versions.items():
        right_flow.append(Paragraph(f"{name}: {ver}", hdr_block))
    right_flow.append(Paragraph(right_sep, hdr_sep))

    header_tbl = Table([[left_flow, right_flow]], colWidths=col_widths_hdr, hAlign="LEFT")
    header_tbl.setStyle(TableStyle([
        ("VALIGN", (0,0), (-1,-1), "TOP"),
        ("ALIGN", (0,0), (-1,-1), "LEFT"),
        ("LEFTPADDING", (0,0), (-1,-1), 0),
        ("RIGHTPADDING", (0,0), (-1,-1), 0),
        ("TOPPADDING", (0,0), (-1,-1), 0),
        ("BOTTOMPADDING", (0,0), (-1,-1), 0),
    ]))
    story.append(header_tbl)
    story.append(Spacer(1, 8))

    # ---------- Render each section ----------
    from bactipipe.relatedness_columns import reference_column_groups, report_sections
    from reportlab.platypus import PageBreak
    from xml.sax.saxutils import escape
    for idx, (sec_title, rows, page_break) in enumerate(report_sections(sections)):
        if page_break:
            story.append(PageBreak())
        heading = sec_title or "Strain Relatedness Summary"
        story.append(Paragraph(heading, section_h))

        if not rows or not rows[0]:
            story.append(Paragraph("<i>No data</i>", hdr_block))
            story.append(Spacer(1, 8))
            continue

        # === 1️⃣ Define column positions here ===
        # === 1️⃣ Define column positions here (alias-aware) ===
        # Extract raw header text (Paragraph or str) and a lowercased copy for matching
        raw_hdr = []
        for h in rows[0]:
            if hasattr(h, "getPlainText"):
                raw_hdr.append(h.getPlainText())
            else:
                raw_hdr.append("" if h is None else str(h))
        hdr_norm = [s.strip().lower() for s in raw_hdr]
        comparison_start, reference_labels = reference_column_groups(raw_hdr)

        def _pos_exact(name: str) -> int:
            """Exact (case-insensitive) lookup against the original header row."""
            key = name.strip().lower()
            try:
                return hdr_norm.index(key)
            except ValueError:
                return -1

        # Tiny alias resolver: try a few label variants without touching the TSV
        def _pos_any(*labels: str) -> int:
            for lab in labels:
                i = _pos_exact(lab)
                if i >= 0:
                    return i
            return -1

        # Column indices (use aliases where headers may vary)
        specimen_idx = _pos_any("specimen")

        size_idx     = _pos_any("size (mb)", "size", "size mb", "genome size", "genome size (mb)")
        gc_idx       = _pos_any("gc (%)", "gc", "gc percent")
        mlst_idx     = _pos_any("mlst", "st")
        cgmlst_idx   = _pos_any("cgmlst", "cgst", "cg-mlst", "cg mlst")
        snp_idx      = _pos_any("snp distance", "snp_dist", "snp distance (ska)", "ska snps vs reference")
        ani_idx      = _pos_any("ani (%)", "ani", "ani percent")

        numeric_cols = {i for i in (size_idx, gc_idx, mlst_idx, cgmlst_idx, snp_idx, ani_idx) if i >= 0}
        if reference_labels:
            numeric_cols.update(range(comparison_start, len(raw_hdr)))


        # === 2️⃣ Then build body_rows (uses specimen_idx & numeric_cols) ===
        wrap_cols = {specimen_idx} if specimen_idx >= 0 else set()
        body_rows: List[List[object]] = []
        for row in rows[1:]:
            formatted_row = []
            for j, val in enumerate(row):
                s = "" if val is None else str(val)
                if j in wrap_cols:
                    formatted_row.append(Paragraph(escape(s), cell_p))
                elif j in numeric_cols:
                    formatted_row.append(Paragraph(f"<para align='center'>{escape(s)}</para>", cell_p))
                else:
                    formatted_row.append(Paragraph(escape(s), cell_p))
            body_rows.append(formatted_row)

        # ----- Canonicalize header *display* labels (PDF only) -----
        canonical_map = {
            "sample": "Sample ID",
            "sample id": "Sample ID",
            "isolate": "Isolate ID",
            "isolate id": "Isolate ID",
            "specimen": "Specimen",
            "size (mb)": "Size (Mb)",
            "size": "Size (Mb)",
            "gc (%)": "GC (%)",
            "gc": "GC (%)",
            "serotype": "Serotype",
            "mlst": "MLST",
            "st": "MLST",
            "cgst": "cgMLST",
            "cgmlst": "cgMLST",
            "snp distance": "SNP<br/>distance",   # two-line header (centered if numeric col)
            "ani (%)": "ANI (%)",
            "ani": "ANI (%)",
        }

        # Raw header texts (strip any Paragraph)
        raw_headers = []
        for h in rows[0]:
            raw_headers.append(h.getPlainText() if isinstance(h, Paragraph) else ("" if h is None else str(h)))

        # Build header cells with per-column alignment
        header_cells = []
        for j, h in enumerate(raw_headers):
            if reference_labels and j >= comparison_start:
                h = 'SNP distance' if (j - comparison_start) % 2 == 0 else 'ANI (%)'
            key = h.strip().lower()
            display = canonical_map.get(key, escape(h))

            # Center only for numeric columns; else left align
            hdr_style = header_cell_center if j in numeric_cols else header_cell_left
            header_cells.append(Paragraph(display, hdr_style))

        # Combine header + body
        table_data: List[List[object]] = [header_cells] + body_rows

        # Column widths: compute and fit
        col_widths_tbl = _column_widths(table_data, doc.width)
        if reference_labels:
            comparison_width = sum(col_widths_tbl[comparison_start:]) / (2 * len(reference_labels))
            col_widths_tbl[comparison_start:] = [comparison_width] * (2 * len(reference_labels))
        spans = []
        header_count = 1
        if reference_labels:
            groups = [''] * len(header_cells)
            groups[0] = Paragraph('Sample information', header_cell_left)
            spans.append(('SPAN', (0, 0), (comparison_start-1, 0)))
            for i, label in enumerate(reference_labels):
                col = comparison_start + 2*i
                groups[col] = Paragraph(escape(label), header_cell_center)
                spans.append(('SPAN', (col, 0), (col+1, 0)))
            table_data.insert(0, groups)
            header_count = 2
        tbl = Table(table_data, colWidths=col_widths_tbl, repeatRows=header_count, hAlign="LEFT")
        tbl.setStyle(TableStyle(spans + [
            ('BACKGROUND', (0, 0), (-1, header_count-1), colors.lightgrey),
        ]))
        tbl.setStyle(TableStyle([
            # Body text
            ("FONT", (0, 1), (-1, -1), "Helvetica", 10),

            # Header row (more readable font for small c/g)
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, 0), 10),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),

            # Background and grid
            ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),

            # Alignment and padding
            ("ALIGN", (0, 0), (-1, -1), "LEFT"),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("TOPPADDING", (0, 0), (-1, -1), 3),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ]))
        story.append(tbl)
        story.append(Spacer(1, 10))

    doc.build(story, onFirstPage=_on_page, onLaterPages=_on_page)



def _read_gc_map_from_assembly(assembly_summary_tsv: Path) -> Dict[str, str]:
    """
    Reads the per-sample GC from the assembly summary if present.
    Returns {sample -> gc_string}. If no GC column, returns {}.
    """
    m: Dict[str, str] = {}
    if not assembly_summary_tsv or not assembly_summary_tsv.exists():
        return m
    # Support either headered or simple columns
    rows = [r for r in csv.reader(assembly_summary_tsv.open("r"), delimiter="\t")]
    if not rows:
        return m
    header = [h.strip().lower() for h in rows[0]]
    if "gc" in header:
        i_sample = header.index("sample")
        i_gc     = header.index("gc")
        for r in rows[1:]:
            if len(r) > max(i_sample, i_gc):
                s = r[i_sample].strip()
                g = r[i_gc].strip()
                if s:
                    m[s] = g or "NA"
    else:
        # no GC column -> nothing to load
        return {}
    return m

def compute_sample_gc_percent(fasta: Path) -> str:
    """
    Streaming GC% across all contigs. Returns formatted string like '52.204'.
    """
    if not fasta or not fasta.exists():
        return ""
    total_len = 0
    gc_count  = 0
    with fasta.open() as fh:
        for line in fh:
            if not line or line.startswith(">"):
                continue
            s = line.strip().upper()
            if not s:
                continue
            total_len += len(s)
            # count G/C quickly
            gc_count += s.count("G") + s.count("C")
    if total_len == 0:
        return ""
    pct = 100.0 * gc_count / total_len
    return f"{pct:.2f}"

# --- matrix builders----------------------------

def _read_tsv_to_map(path: Path, key_col: int = 0, val_col: int = 1) -> Dict[str, str]:
    m: Dict[str, str] = {}
    if not path or not path.exists():
        return m
    with path.open("r") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            if len(row) <= max(key_col, val_col):
                continue
            k = row[key_col].strip()
            v = row[val_col].strip()
            if k and k.lower() != "sample":
                m[k] = v
    return m

def triangular_matrix_from_pairs(
    pairs_tsv: Path,
    sample_order: List[str],
    *,
    sample1_keys: Iterable[str] = ("Sample1", "sample1", "Ref", "ref", "File1"),
    sample2_keys: Iterable[str] = ("Sample2", "sample2", "Query", "query", "File2"),
    value_keys:    Iterable[str] = ("ANI", "ani", "Distance", "distance"),
    diagonal_value: str = "0",
    placeholder:    str = "-",
    value_format: Optional[str] = None,   # e.g. "{:.2f}" or "{:,}"
) -> List[List[str]]:
    """
    Build a strictly LOWER-TRIANGULAR matrix:
      - Diagonal = diagonal_value (e.g., "0")
      - Lower triangle = values from pairs_tsv
      - Upper triangle = placeholder (e.g., "-")
    Pairs are treated as unordered (A,B) == (B,A); we prefer numeric if both exist.
    """
    if not pairs_tsv or not pairs_tsv.exists():
        return []

    rows = [r for r in csv.reader(pairs_tsv.open("r"), delimiter="\t") if any(c.strip() for c in r)]
    if not rows:
        return []

    header = [h.strip() for h in rows[0]]
    hmap = {h.lower(): i for i, h in enumerate(header)}

    def _first_index(keys: Iterable[str]) -> Optional[int]:
        for k in keys:
            idx = hmap.get(k.lower())
            if idx is not None:
                return idx
        return None

    i_s1 = _first_index(sample1_keys)
    i_s2 = _first_index(sample2_keys)
    i_v  = _first_index(value_keys)
    if i_s1 is None or i_s2 is None or i_v is None:
        return []

    # Map unordered pair -> value
    def _coerce_num(s: str) -> Optional[float]:
        try:
            return float(s)
        except Exception:
            return None

    d: Dict[frozenset, str] = {}
    for r in rows[1:]:
        if len(r) <= max(i_s1, i_s2, i_v):
            continue
        a, b, v = r[i_s1].strip(), r[i_s2].strip(), r[i_v].strip()
        if not a or not b:
            continue
        key = frozenset((a, b))
        # Keep the better (numeric) value if we see the pair twice
        if key in d:
            old_num = _coerce_num(d[key])
            new_num = _coerce_num(v)
            if old_num is None and new_num is not None:
                d[key] = v
        else:
            d[key] = v

    # Now build matrix
    out: List[List[str]] = []
    out.append(["sample"] + sample_order)
    for i, si in enumerate(sample_order):
        row = [si]
        for j, sj in enumerate(sample_order):
            if i == j:
                row.append(diagonal_value)
            elif i > j:
                key = frozenset((si, sj))
                val = d.get(key, placeholder)
                if value_format and val not in {placeholder, ""}:
                    try:
                        val = value_format.format(float(val))
                    except Exception:
                        pass
                row.append(val)
            else:
                row.append(placeholder)
        out.append(row)
    return out

# --- main writer --------------------------------------

def write_final_summary(
    *,
    accession: str,
    manifest: List[SampleRecord],
    paths: Dict[str, Path],
    logger,
    reference_provided: Optional[bool] = None,
    has_reference: Optional[bool] = None,   # alias
    reference_comparisons=None,
) -> Path:
    """
    Compose the consolidated final TSV at:
      paths['outdir'] / f"{accession}_relate.tsv"

    Summary columns (with reference):
      Sample | Isolate | Specimen | Size (Mb) | GC (%) | Serotype | ST | cgST | SNP distance | ANI (%)

    When NO reference:
      Drop ANI and SNP distance from the summary table, and append LOWER-triangular:
        - '## ANI matrix (skani)'
        - '## SKA SNPs matrix'
    """
    if has_reference is None:
        has_reference = bool(reference_provided)
    comparisons = reference_comparisons
    if comparisons is None:
        comparisons = ([dict(sample_id=manifest[0].sample, ani_tsv=paths.get('ani_tsv'),
                             ska_vs_ref=paths.get('ska_vs_ref'))] if has_reference and manifest else [])
    has_reference = bool(comparisons)

    outdir = paths["outdir"]
    outdir.mkdir(parents=True, exist_ok=True)
    final_path = outdir / f"{accession}_relate.tsv"

    # --- local helpers ---
    def _read_size_map_from_assembly(asm_tsv: Optional[Path]) -> Dict[str, str]:
        """Return {sample: 'X.XXX'} from 'size_mb' or compute from 'sum_len'."""
        out: Dict[str, str] = {}
        if not asm_tsv or not asm_tsv.exists():
            return out
        try:
            with asm_tsv.open("r") as fh:
                reader = csv.reader(fh, delimiter="\t")
                header = next(reader, None)
                if not header:
                    return out
                h = [c.strip().lower() for c in header]
                i_s = h.index("sample") if "sample" in h else None
                i_mb = h.index("size_mb") if "size_mb" in h else None
                i_bp = h.index("sum_len") if "sum_len" in h else None
                for row in reader:
                    if not row or i_s is None or i_s >= len(row):
                        continue
                    s = row[i_s].strip()
                    if i_mb is not None and i_mb < len(row) and row[i_mb].strip():
                        out[s] = row[i_mb].strip()
                    elif i_bp is not None and i_bp < len(row):
                        try:
                            bp = float(row[i_bp])
                            out[s] = f"{bp/1_000_000.0:.2f}"
                        except Exception:
                            pass
        except Exception as e:
            logger.warning(f"Failed reading assembly summary for size_mb: {e}")
        return out

    def _fmt_int_commas_or_dash(val: Optional[str]) -> str:
        """SNPs: integer with commas; '-' passthrough; else 'NA'."""
        if val is None:
            return "NA"
        v = str(val).strip()
        if v in {"", "NA", "na", "NaN"}:
            return "NA"
        if v == "-":
            return "-"
        try:
            f = float(v.replace(",", ""))
            return f"{int(round(f)):,}"
        except Exception:
            return "NA"

    def _fmt_pct_or_dash(val: Optional[str], nd=2) -> str:
        """Percent: '99.863'→'99.86'; '-' passthrough; else 'NA'."""
        if val is None:
            return "NA"
        v = str(val).strip()
        if v == "-":
            return "-"
        try:
            f = float(v)
            return f"{f:.{nd}f}"
        except Exception:
            return "NA"

    # --- maps from upstream steps ---
    sero_map = _read_tsv_to_map(paths.get("serotypes_tsv")) if paths.get("serotypes_tsv") else {}
    st_map   = _read_tsv_to_map(paths.get("mlst_tsv"))      if paths.get("mlst_tsv")      else {}
    cgst_map = _read_tsv_to_map(paths.get("cgmlst_tsv"))    if paths.get("cgmlst_tsv")    else {}
    comparison_maps = []
    for comparison in comparisons:
        def metric_map(key):
            path = comparison.get(key)
            return _read_tsv_to_map(path) if path and path.exists() else {}
        comparison_maps.append((comparison['sample_id'], metric_map('ska_vs_ref'), metric_map('ani_tsv')))

    # GC (%) and Size (Mb) from assembly summary
    gc_map   = _read_gc_map_from_assembly(paths.get("assembly_summary_tsv")) if paths.get("assembly_summary_tsv") else {}
    size_map = _read_size_map_from_assembly(paths.get("assembly_summary_tsv"))

    # --- write final TSV ---
    with final_path.open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")

        # Title row
        w.writerow(["## Strain Relatedness Summary"])

        header = ["Sample", "Isolate", "Specimen", "Size (Mb)", "GC (%)", "Serotype", "ST", "cgST"]
        for ref_id, _, _ in comparison_maps:
            prefix = f'{ref_id}: ' if len(comparisons) > 1 else ''
            header += [prefix + 'SNP distance', prefix + 'ANI (%)']
        w.writerow(header)

        ref_sample = manifest[0].sample if (has_reference and manifest) else None

        for rec in manifest:
            if not rec.exists:
                continue
            s   = rec.sample
            iso = rec.isolate
            spe = rec.specimen

            size_mb = size_map.get(s, "NA")
            gc_pct  = gc_map.get(s, "NA")
            sero    = sero_map.get(s, "NA")
            st      = st_map.get(s, "NA")
            cgst    = cgst_map.get(s, "NA")

            row = [s, iso, spe, size_mb, gc_pct, sero, st, cgst]

            for ref_id, ska_vs_ref_map, ani_map in comparison_maps:
                is_ref = (s == ref_id)
                ani_raw = "-" if is_ref else ani_map.get(s, "NA")
                snp_raw = "-" if is_ref else ska_vs_ref_map.get(s, "NA")
                ani_fmt = _fmt_pct_or_dash(ani_raw, nd=2)
                snp_fmt = _fmt_int_commas_or_dash(snp_raw)
                row += [snp_fmt, ani_fmt]

            w.writerow(row)

        w.writerow([])  # spacer

        # Matrices only when NO reference
        if not has_reference:
            sample_order = [r.sample for r in manifest if r.exists]

            # ANI triangular (skani)
            ani_pairs = paths.get("skani_pairs_samples") or paths.get("skani_pairs")
            ani_mat = triangular_matrix_from_pairs(
                ani_pairs, sample_order,
                sample1_keys=("Sample1","sample1","Ref","ref"),
                sample2_keys=("Sample2","sample2","Query","query"),
                value_keys=("ANI","ani"),
                diagonal_value="0",
                placeholder="-",
                value_format="{:.2f}",
            )
            if ani_mat:
                w.writerow(["## ANI matrix (skani)"])
                for r in ani_mat:
                    w.writerow(r)
                w.writerow([])

            # SKA SNPs triangular (Distance), with thousands-separators
            ska_mat = triangular_matrix_from_pairs(
                paths.get("ska_distances"), sample_order,
                sample1_keys=("Sample1","sample1","Ref","ref"),
                sample2_keys=("Sample2","sample2","Query","query"),
                value_keys=("Distance","distance"),
                diagonal_value="0",
                placeholder="-",
                value_format="{:,}",
            )
            if ska_mat:
                w.writerow(["## SKA SNPs matrix"])
                for r in ska_mat:
                    w.writerow(r)

    logger.info(f"Wrote final summary: {final_path}")
    return final_path


def run_skani_triangle_and_write(
    samples: List[Tuple[str, Path]],   # [(sample_name, fasta_path)]
    skani_root: Path,
    out_prefix: str,
    logger,
    env_name: Optional[str] = None,
) -> Dict[str, Path]:
    """
    Calls:  skani triangle -o <triangle.tsv> <FASTA...>
    Parses its upper-triangle output into:
      - <prefix>.triangle.tsv        (raw triangle file from skani)
      - <prefix>.pairs.tsv           (by file paths)
      - <prefix>.pairs_samples.tsv   (by sample names)
    Returns a dict with those paths (only the ones that exist).
    """
    skani_root.mkdir(parents=True, exist_ok=True)
    # Use absolute paths so triangle output paths match exactly
    fasta_paths = [str(p.resolve()) for _, p in samples]
    tri_path   = skani_root / f"{out_prefix}.triangle.tsv"
    stdout_path = skani_root / f"{out_prefix}.triangle.stdout.txt"

    cmd = ["skani", "triangle", "-o", str(tri_path)] + fasta_paths
    rc, out = _run(cmd, cwd=None, log=logger, env_name=env_name)
    stdout_path.write_text(out)

    if rc != 0 or not tri_path.exists():
        logger.error("skani triangle FAILED (no triangle file). See %s", stdout_path)
        return {}

    # Parse triangle format
    lines = [ln.rstrip("\n") for ln in tri_path.open("r") if ln.strip()]
    try:
        n = int(lines[0].split()[0])
    except Exception:
        logger.error("skani triangle parse error: first line does not contain N. See %s", tri_path)
        return {"triangle": tri_path}

    # Map triangle file paths -> sample names
    name_by_file: Dict[str, str] = {str(p.resolve()): s for s, p in samples}

    files_in_order: List[str] = []
    pairs_file = skani_root / f"{out_prefix}.pairs.tsv"
    pairs_samples_file = skani_root / f"{out_prefix}.pairs_samples.tsv"

    with pairs_file.open("w", newline="") as pf, pairs_samples_file.open("w", newline="") as psf:
        pw = csv.writer(pf, delimiter="\t")
        sw = csv.writer(psf, delimiter="\t")
        pw.writerow(["File1", "File2", "ANI"])
        sw.writerow(["Sample1", "Sample2", "ANI"])

        rows = lines[1:1+n]
        for i, row in enumerate(rows):
            cols = row.split("\t")
            fi = cols[0].strip()
            files_in_order.append(fi)
            vals = [c.strip() for c in cols[1:] if c.strip() != ""]
            # vals are ANI between current file (i) and previous files j < i
            for j, v in enumerate(vals):
                f1 = files_in_order[j]
                f2 = fi
                try:
                    ani = f"{float(v):.2f}"
                except Exception:
                    ani = "NA"
                pw.writerow([f1, f2, ani])
                s1 = name_by_file.get(f1, Path(f1).stem)
                s2 = name_by_file.get(f2, Path(f2).stem)
                sw.writerow([s1, s2, ani])

    return {
        "triangle": tri_path,
        "pairs": pairs_file,
        "pairs_samples": pairs_samples_file,
    }
