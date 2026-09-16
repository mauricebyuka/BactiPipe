#!/usr/bin/env python3
"""Build a reviewed-name, chromosome-only genome-size catalog from NCBI.

The script never changes the controlled organism vocabulary. It queries RefSeq
assemblies for each taxid, retrieves the assembly sequence reports, sums only
sequences explicitly assigned to chromosomes, and publishes an auditable
candidate TSV. Network stages are checkpointed and safe to resume.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import statistics
import subprocess
import time


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def status(path: Path, stage: str, state: str, detail: str = "") -> None:
    timestamp = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    atomic_text(path, f"{timestamp}\t{stage}\t{state}\t{detail}\n")


def controlled_rows(path: Path) -> list[dict[str, object]]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                raise ValueError(f"Invalid controlled row {line_number}.")
            rows.append(
                {
                    "taxid": parts[0].strip(),
                    "legacy_size_bp": int(parts[1]),
                    "controlled_accession": parts[2].strip(),
                    "canonical_name": parts[-1].strip(),
                }
            )
    return rows


def run_datasets(arguments: list[str], attempts: int = 4) -> str:
    delay = 2.0
    for attempt in range(1, attempts + 1):
        result = subprocess.run(
            ["datasets", *arguments],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if result.returncode == 0:
            return result.stdout
        if attempt == attempts:
            raise RuntimeError(result.stderr.strip() or f"datasets exit {result.returncode}")
        time.sleep(delay)
        delay *= 2
    raise AssertionError("unreachable")


def json_lines(text: str) -> list[dict[str, object]]:
    return [json.loads(line) for line in text.splitlines() if line.strip().startswith("{")]


def chunks(items: list[str], size: int) -> list[list[str]]:
    return [items[index : index + size] for index in range(0, len(items), size)]


def percentile(values: list[int], fraction: float) -> int:
    if len(values) == 1:
        return values[0]
    position = (len(values) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(values) - 1)
    weight = position - lower
    return int(round(values[lower] * (1 - weight) + values[upper] * weight))


def summarize_catalog(
    controlled: list[dict[str, object]],
    accession_taxids: dict[str, set[str]],
    sequence_reports: list[dict[str, object]],
) -> list[dict[str, object]]:
    chromosome_by_assembly: dict[str, int] = {}
    for report in sequence_reports:
        if str(report.get("assigned_molecule_location_type", "")).casefold() != "chromosome":
            continue
        accession = str(report.get("assembly_accession", "")).strip()
        length = int(report.get("length", 0) or 0)
        if accession and length > 0:
            chromosome_by_assembly[accession] = (
                chromosome_by_assembly.get(accession, 0) + length
            )

    sizes_by_taxid: dict[str, list[tuple[str, int]]] = {}
    for accession, taxids in accession_taxids.items():
        size = chromosome_by_assembly.get(accession)
        if not size:
            continue
        for taxid in taxids:
            sizes_by_taxid.setdefault(taxid, []).append((accession, size))

    result = []
    for row in controlled:
        taxid = str(row["taxid"])
        observations = sorted(sizes_by_taxid.get(taxid, []))
        if observations:
            values = sorted(size for _, size in observations)
            result.append(
                {
                    **row,
                    "expected_size_bp": int(round(statistics.median(values))),
                    "lower_size_bp": percentile(values, 0.10),
                    "upper_size_bp": percentile(values, 0.90),
                    "assembly_count": len(values),
                    "method": "median_refseq_chromosome_length",
                    "source_accessions": ",".join(
                        accession for accession, _ in observations
                    ),
                }
            )
        else:
            result.append(
                {
                    **row,
                    "expected_size_bp": row["legacy_size_bp"],
                    "lower_size_bp": "N/A",
                    "upper_size_bp": "N/A",
                    "assembly_count": 0,
                    "method": "controlled_reference_fallback_no_chromosome_report",
                    "source_accessions": row["controlled_accession"],
                }
            )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--controlled-list", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--checkpoint-dir", type=Path, required=True)
    parser.add_argument("--status-file", type=Path, required=True)
    parser.add_argument("--assemblies-per-taxon", type=int, default=25)
    parser.add_argument("--sequence-batch-size", type=int, default=100)
    parser.add_argument("--request-delay", type=float, default=0.35)
    args = parser.parse_args()

    args.checkpoint_dir.mkdir(parents=True, exist_ok=True)
    assembly_dir = args.checkpoint_dir / "assemblies"
    sequence_dir = args.checkpoint_dir / "sequences"
    assembly_dir.mkdir(exist_ok=True)
    sequence_dir.mkdir(exist_ok=True)
    failures_path = args.checkpoint_dir / "query_failures.tsv"
    controlled = controlled_rows(args.controlled_list)
    status(args.status_file, "assemblies", "running", f"taxa={len(controlled)}")

    accession_taxids: dict[str, set[str]] = {}
    failures: list[tuple[str, str, str]] = []
    for index, row in enumerate(controlled, start=1):
        taxid = str(row["taxid"])
        checkpoint = assembly_dir / f"{taxid}.jsonl"
        if not checkpoint.exists():
            try:
                output = run_datasets(
                    [
                        "summary",
                        "genome",
                        "taxon",
                        taxid,
                        "--assembly-level",
                        "complete,chromosome",
                        "--assembly-source",
                        "RefSeq",
                        "--exclude-atypical",
                        "--exclude-multi-isolate",
                        "--mag",
                        "exclude",
                        "--limit",
                        str(args.assemblies_per_taxon),
                        "--as-json-lines",
                    ]
                )
                atomic_text(checkpoint, output)
            except RuntimeError as exc:
                failures.append((taxid, "assembly_summary", str(exc).replace("\t", " ")))
                continue
            time.sleep(max(0.0, args.request_delay))
        for report in json_lines(checkpoint.read_text(encoding="utf-8")):
            accession = str(report.get("accession", "")).strip()
            if accession:
                accession_taxids.setdefault(accession, set()).add(taxid)
        if index % 100 == 0:
            status(args.status_file, "assemblies", "running", f"{index}/{len(controlled)}")

    accessions = sorted(accession_taxids)
    status(args.status_file, "sequences", "running", f"assemblies={len(accessions)}")
    sequence_reports: list[dict[str, object]] = []
    for batch in chunks(accessions, args.sequence_batch_size):
        digest = sha256("\n".join(batch).encode("utf-8")).hexdigest()[:16]
        checkpoint = sequence_dir / f"{digest}.jsonl"
        input_path = sequence_dir / f"{digest}.accessions.txt"
        if not checkpoint.exists():
            atomic_text(input_path, "\n".join(batch) + "\n")
            try:
                output = run_datasets(
                    [
                        "summary",
                        "genome",
                        "accession",
                        "--inputfile",
                        str(input_path),
                        "--report",
                        "sequence",
                        "--as-json-lines",
                    ]
                )
                atomic_text(checkpoint, output)
            except RuntimeError as exc:
                failures.append((digest, "sequence_report", str(exc).replace("\t", " ")))
                continue
            time.sleep(max(0.0, args.request_delay))
        sequence_reports.extend(json_lines(checkpoint.read_text(encoding="utf-8")))

    rows = summarize_catalog(controlled, accession_taxids, sequence_reports)
    fields = [
        "taxid",
        "canonical_name",
        "expected_size_bp",
        "lower_size_bp",
        "upper_size_bp",
        "assembly_count",
        "method",
        "source_accessions",
        "controlled_accession",
        "legacy_size_bp",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, dialect="excel-tab")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(args.output)

    with failures_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, dialect="excel-tab")
        writer.writerow(["identifier", "stage", "error"])
        writer.writerows(failures)
    status(
        args.status_file,
        "catalog",
        "success",
        f"rows={len(rows)}; chromosome_sizes={sum(int(row['assembly_count']) > 0 for row in rows)}; failures={len(failures)}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
