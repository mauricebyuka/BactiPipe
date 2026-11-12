"""Pure execution-planning helpers shared by BactiPipe QC runners."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from bactipipe.genome_sizes import GenomeSizeResolution


@dataclass(frozen=True)
class NanoporeQCJob:
    sample: str
    organism_entered: str
    organism_canonical: str
    barcode: str
    raw_folder: str
    output_dir: str
    output_fastq: str
    genome_size_bp: int | None
    min_avg_quality: float
    min_read_length: int
    min_coverage: float
    desired_coverage: float
    cpus: int
    use_s3: bool = False
    bucket: str | None = None


@dataclass(frozen=True)
class NanoporeAssemblyJob:
    sample: str
    fastq: str
    genome_size_bp: int | None
    assembly_dir: str
    assembler: str
    cpus: int
    logfile: str
    min_quality: float = 15
    min_length: int = 1000


def bounded_worker_count(
    *,
    sample_count: int,
    total_cpus: int,
    cpus_per_sample: int,
    memory_worker_limit: int,
) -> int:
    """Bound workers by samples, CPU allocation, and the memory estimate."""

    if sample_count <= 0:
        return 0
    total_cpus = max(1, int(total_cpus))
    cpus_per_sample = max(1, int(cpus_per_sample))
    return max(
        1,
        min(
            int(sample_count),
            total_cpus // cpus_per_sample,
            max(1, int(memory_worker_limit)),
        ),
    )


def flye_size_arguments(genome_size_bp: int | None) -> list[str]:
    """Only known sizes enable Flye's size-dependent initial downsampling."""

    if genome_size_bp is None:
        return []
    return ["--asm-coverage", "50", "--g", str(genome_size_bp)]


def nanopore_qc_failure_reason(
    metrics_text: str,
    *,
    genome_size_bp: int | None,
    minimum_quality: float,
    minimum_coverage: float,
) -> str | None:
    """Reject unusable reads, not low predicted depth.

    Genome size and minimum coverage remain accepted for caller compatibility.
    Depth is a reporting threshold assessed against the actual assembly, not a
    pre-assembly execution gate based on the expected organism's genome size.
    """

    def metric(label: str) -> float | None:
        match = re.search(
            rf"^{re.escape(label)}\s*:\s*([0-9][0-9,.]*)\s*$",
            metrics_text,
            flags=re.MULTILINE,
        )
        if not match:
            return None
        return float(match.group(1).replace(",", ""))

    total_bases = metric("Trimmed Total Bases")
    mean_quality = metric("Trimmed Mean Read Quality")
    if total_bases is None or mean_quality is None:
        return "Read quality filtering did not produce complete QC outputs."
    if total_bases <= 0:
        return "Read quality filtering retained no sequence bases."
    if mean_quality < minimum_quality:
        return (
            f"Trimmed mean read quality {mean_quality:.2f} is below the required "
            f"{minimum_quality:g}."
        )
    return None


def illumina_minimum_coverage(mean_quality: float, canonical_organism: str | None) -> int:
    """Return the established BactiPipe Illumina coverage threshold."""

    if mean_quality >= 30:
        minimum = 30
    elif mean_quality >= 29:
        minimum = 40
    else:
        minimum = 50
    if canonical_organism == "Escherichia coli":
        minimum += 10
    return minimum


def overall_quality_verdict(
    *,
    quality_verdict: str,
    coverage_verdict: str,
    taxonomy_verdict: str,
    completeness_verdict: str,
    contamination_verdict: str,
) -> str:
    """Require measured QC metrics while allowing no expected taxonomy."""

    from bactipipe.qc_status import quality_status
    quality = quality_status([quality_verdict, coverage_verdict, completeness_verdict, contamination_verdict])
    if quality == "Fail":
        return "Fail"
    return "Review" if quality != "Pass" or taxonomy_verdict in {"Fail", "Review"} else "Pass"


def write_execution_plan(
    path: str | Path,
    rows: Iterable[tuple[str, str, GenomeSizeResolution]],
    *,
    platform: str,
    min_quality: float,
    min_read_length: int | None,
    min_coverage: float | str,
    desired_coverage: float | None = None,
    catalog_checksum: str = "N/A",
    assembly_uses_genome_size: bool = True,
) -> None:
    """Persist the exact per-sample size decisions used by the run."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    fieldnames = [
        "sample",
        "entered_organism",
        "canonical_organism",
        "resolution_status",
        "taxid",
        "genome_size_bp",
        "size_source",
        "reference_accession",
        "catalog_version",
        "catalog_sha256",
        "platform",
        "minimum_read_quality",
        "minimum_read_length",
        "minimum_coverage_depth",
        "qc_target_bases",
        "assembly_size_mode",
    ]
    with tmp_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, dialect="excel-tab")
        writer.writeheader()
        for sample, entered, resolution in rows:
            genome_size = resolution.genome_size_bp
            writer.writerow(
                {
                    "sample": sample,
                    "entered_organism": entered,
                    "canonical_organism": resolution.canonical_name or "N/A",
                    "resolution_status": resolution.status,
                    "taxid": resolution.taxid or "N/A",
                    "genome_size_bp": genome_size if genome_size is not None else "N/A",
                    "size_source": resolution.source or "N/A",
                    "reference_accession": resolution.reference_accession or "N/A",
                    "catalog_version": resolution.catalog_version,
                    "catalog_sha256": catalog_checksum,
                    "platform": platform,
                    "minimum_read_quality": min_quality,
                    "minimum_read_length": (
                        min_read_length if min_read_length is not None else "N/A"
                    ),
                    "minimum_coverage_depth": min_coverage,
                    "qc_target_bases": (
                        int(genome_size * desired_coverage)
                        if genome_size is not None and desired_coverage is not None
                        else "N/A"
                    ),
                    "assembly_size_mode": (
                        "known_genome_size"
                        if assembly_uses_genome_size and genome_size is not None
                        else "size_independent"
                    ),
                }
            )
    tmp_path.replace(path)
