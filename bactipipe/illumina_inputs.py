"""Resolve physical Illumina inputs independently of the output run name."""
import json
from pathlib import Path


def stage_s3_reads(*, bucket, prefix, run_name, input_run_name, target_root, download):
    source_name = input_run_name or run_name
    if source_name in {".", ".."} or "/" in source_name or "\\" in source_name:
        raise ValueError("The input run name must be a single directory name.")
    source_prefix = "/".join(
        part.strip("/") for part in (prefix, "DATA", "illuminaData", source_name)
        if part and part.strip("/")
    )
    staged = Path(target_root) / source_name
    marker = staged / ".bactipipe-download-complete.json"
    identity = {"bucket": bucket, "prefix": source_prefix}
    try:
        reusable = json.loads(marker.read_text()) == identity
    except (OSError, ValueError):
        reusable = False
    if not reusable:
        download(bucket=bucket, prefix=source_prefix, target_root=str(target_root))
    if not any(
        p.is_file() and p.stat().st_size > 0
        for p in staged.rglob("*")
        if p.name.lower().endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))
    ):
        raise ValueError(f"No FASTQ files found at s3://{bucket}/{source_prefix}/")
    if not reusable:
        temporary = marker.with_suffix(".tmp")
        temporary.write_text(json.dumps(identity) + "\n")
        temporary.replace(marker)
    return str(staged)
