"""Optional Kraken2 candidate and mixed-signal evidence collection."""

from __future__ import annotations

import csv
import os
import subprocess
from pathlib import Path
from typing import Mapping, Sequence

from bactipipe.species_identification import TaxonomicSignal


def run_kraken2(
    *,
    reads: Sequence[str | Path],
    database: str | Path,
    output_directory: str | Path,
    sample_id: str,
    threads: int = 1,
    confidence: float = 0.05,
    executable: str = "kraken2",
) -> Path:
    """Run Kraken2 with distinct-minimizer reporting.

    The output nominates candidate taxa.  It does not itself constitute a
    positive species identification or prove a mixed culture.
    """

    if len(reads) not in {1, 2}:
        raise ValueError("Kraken2 evidence requires one read file or an Illumina pair")
    read_paths = [Path(path) for path in reads]
    for path in read_paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Read file is missing or empty: {path}")
    database_path = Path(database)
    if not database_path.exists():
        raise ValueError(f"Kraken2 database does not exist: {database_path}")
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    report = output / "kraken2.report.tsv"
    classifications = output / "kraken2.classifications.tsv"
    temporary_report = report.with_suffix(report.suffix + ".tmp")
    temporary_classifications = classifications.with_suffix(
        classifications.suffix + ".tmp"
    )
    command = [
        executable,
        "--db",
        str(database_path),
        "--threads",
        str(max(1, int(threads))),
        "--confidence",
        str(confidence),
        "--report-minimizer-data",
        "--report",
        str(temporary_report),
        "--output",
        str(temporary_classifications),
    ]
    if len(read_paths) == 2:
        command.append("--paired")
    command.extend(str(path) for path in read_paths)
    try:
        completed = subprocess.run(command, text=True, capture_output=True)
    except OSError as exc:
        temporary_report.unlink(missing_ok=True)
        temporary_classifications.unlink(missing_ok=True)
        raise RuntimeError(f"Kraken2 could not be started for {sample_id}: {exc}") from exc
    (output / "kraken2.stderr.log").write_text(completed.stderr or "")
    if completed.returncode != 0:
        temporary_report.unlink(missing_ok=True)
        temporary_classifications.unlink(missing_ok=True)
        raise RuntimeError(
            f"Kraken2 failed for {sample_id} with exit {completed.returncode}; "
            f"see {output / 'kraken2.stderr.log'}"
        )
    if not temporary_report.is_file():
        raise RuntimeError("Kraken2 completed without producing its report")
    os.replace(temporary_report, report)
    if temporary_classifications.exists():
        os.replace(temporary_classifications, classifications)
    return report


def parse_kraken2_report(
    path: str | Path,
    *,
    identification_names_by_taxid: Mapping[str, str] | None = None,
) -> tuple[TaxonomicSignal, ...]:
    """Parse candidate-taxon counts from an eight-column minimizer report.

    When a reviewed identification catalog is supplied, its callable taxids
    define the reporting units even when NCBI ranks one of them below species.
    Direct rather than clade counts then prevent a parent species row from
    double-counting a retained clinical subspecies row.
    """

    signals: list[TaxonomicSignal] = []
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 8:
                raise ValueError(
                    f"Kraken2 minimizer report line {line_number} has {len(fields)} "
                    "columns; expected 8"
                )
            percent, clade_reads, direct_reads, _minimizers, distinct, rank, taxid, name = fields
            taxid = taxid.strip()
            if identification_names_by_taxid is None:
                if rank != "S":
                    continue
                selected_reads = int(clade_reads)
                selected_percent = float(percent)
                selected_name = " ".join(name.split())
            else:
                selected_name = identification_names_by_taxid.get(taxid, "")
                if not selected_name:
                    continue
                selected_reads = int(direct_reads)
                if selected_reads <= 0:
                    continue
                selected_percent = (
                    float(percent) * selected_reads / int(clade_reads)
                    if int(clade_reads)
                    else 0.0
                )
            signals.append(
                TaxonomicSignal(
                    species_taxid=taxid,
                    species_name=selected_name,
                    read_count=selected_reads,
                    read_percent=selected_percent,
                    reference_breadth=None,
                    distinct_minimizers=int(distinct),
                )
            )
    signals.sort(key=lambda item: (item.read_count, item.read_percent), reverse=True)
    return tuple(signals)


def write_taxonomic_signals(
    signals: Sequence[TaxonomicSignal], path: str | Path
) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    fields = [
        "species_taxid",
        "species_name",
        "read_count",
        "read_percent",
        "distinct_minimizers",
        "reference_breadth",
        "competitive_mapped_reads",
        "competitive_mapped_percent",
        "nomination_score",
        "nomination_percent",
    ]
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for signal in signals:
            writer.writerow(
                {
                    "species_taxid": signal.species_taxid,
                    "species_name": signal.species_name,
                    "read_count": signal.read_count,
                    "read_percent": f"{signal.read_percent:.4f}",
                    "distinct_minimizers": (
                        signal.distinct_minimizers
                        if signal.distinct_minimizers is not None
                        else "N/A"
                    ),
                    "reference_breadth": (
                        f"{signal.reference_breadth:.4f}"
                        if signal.reference_breadth is not None
                        else "N/A"
                    ),
                    "competitive_mapped_reads": (
                        signal.competitive_mapped_reads
                        if signal.competitive_mapped_reads is not None
                        else "N/A"
                    ),
                    "competitive_mapped_percent": (
                        f"{signal.competitive_mapped_percent:.4f}"
                        if signal.competitive_mapped_percent is not None
                        else "N/A"
                    ),
                    "nomination_score": (
                        f"{signal.nomination_score:.4f}"
                        if signal.nomination_score is not None
                        else "N/A"
                    ),
                    "nomination_percent": (
                        f"{signal.nomination_percent:.4f}"
                        if signal.nomination_percent is not None
                        else "N/A"
                    ),
                }
            )
    os.replace(temporary, destination)
