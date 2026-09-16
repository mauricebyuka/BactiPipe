"""Per-sample failure audit helpers for BactiPipe QC batches."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


FAILURE_COLUMNS = ("sample_id", "stage", "reason")


def record_failure(
    failures: dict[str, tuple[str, str]],
    sample_id: str,
    stage: str,
    reason: object,
) -> None:
    """Keep the first actionable failure for a sample."""

    sample_id = str(sample_id or "").strip()
    if not sample_id or sample_id in failures:
        return
    clean_reason = " ".join(str(reason or "Sample analysis failed.").split())
    failures[sample_id] = (str(stage or "analysis").strip(), clean_reason)


def publish_failures(
    outdir: str | Path,
    temp_summary: str | Path,
    expected_organisms: dict[str, str],
    failures: dict[str, tuple[str, str]],
    *,
    required_depth: int | float,
) -> Path:
    """Write the audit and add one parseable placeholder row per failed sample."""

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    failure_path = outdir / "sample_failures.tsv"
    with failure_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(FAILURE_COLUMNS)
        for sample_id, (stage, reason) in failures.items():
            writer.writerow((sample_id, stage, reason))

    temp_summary = Path(temp_summary)
    existing = set()
    if temp_summary.is_file():
        with temp_summary.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if row:
                    existing.add(row[0])

    with temp_summary.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for sample_id in expected_organisms:
            if sample_id not in failures or sample_id in existing:
                continue
            depth, depth_verdict, threshold = "N/A", "N/A", required_depth
            # If identification failed after successful assembly, preserve measured
            # depth while retaining the sample's independent execution failure.
            if failures[sample_id][0] == "organism_identification":
                try:
                    audit = json.loads((outdir / "qc_out" / sample_id / "depth_metrics.json").read_text())
                    value = float(audit["depth"])
                    if math.isfinite(value) and value >= 0:
                        depth = f"{value:.2f}X"
                        depth_verdict = audit["verdict"]
                        threshold = audit["minimum_depth"]
                except (OSError, ValueError, TypeError, KeyError):
                    pass
            writer.writerow(
                (
                    sample_id,
                    "0",
                    "Fail",
                    expected_organisms[sample_id],
                    "Analysis failed",
                    "0",
                    depth,
                    threshold,
                    depth_verdict,
                    "Fail",
                )
            )
    return failure_path
