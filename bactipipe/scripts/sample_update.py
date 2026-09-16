"""Transactional selected-sample updates for BactiPipe QC commands."""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


def run_qc_sample_update(
    *,
    script_path: str | Path,
    argv: list[str],
    outdir_parent: str | Path,
    run_name: str,
    sample_names: list[str],
    species_database_root: Optional[str] = None,
    species_platform: str = "",
    species_gtdb_metadata: Optional[str] = None,
    species_gtdbtk_conda_env: str = "gtdbtk",
    species_placement_threads: int = 1,
) -> None:
    """Run selected samples separately, merge their outputs, and publish atomically."""

    selected = tuple(dict.fromkeys(name.strip() for name in sample_names if name.strip()))
    if not selected:
        raise ValueError("--update-existing requires at least one --only-sample value.")
    parent = Path(outdir_parent).resolve()
    target = parent / run_name
    _validate_existing_results(target, run_name)
    stage_parent = Path(
        tempfile.mkdtemp(prefix=f".{run_name}.sample-update-", dir=parent)
    )
    fresh_parent = stage_parent / "fresh"
    fresh = fresh_parent / run_name
    merged = stage_parent / "merged"
    fresh_parent.mkdir()
    try:
        command = [
            sys.executable,
            str(Path(script_path).resolve()),
            *_replace_outdir_argument(argv, fresh_parent),
        ]
        completed = subprocess.run(command, check=False)
        if completed.returncode:
            raise RuntimeError(
                "Selected-sample update was not published because reanalysis exited "
                f"with status {completed.returncode}."
            )
        failures = _selected_failures(fresh / "sample_failures.tsv", selected)
        if failures:
            details = "; ".join(
                f"{row['sample_id']}: {row.get('reason') or 'analysis failed'}"
                for row in failures
            )
            raise RuntimeError(
                "Selected-sample update was not published because reanalysis failed: "
                + details
            )

        shutil.copytree(target, merged)
        for sample_name in selected:
            _replace_sample_outputs(fresh, merged, sample_name)
        _merge_qc_summary(target, fresh, merged, run_name, selected)
        _merge_temporary_summary(target, fresh, merged, selected)
        _merge_failure_audit(target, fresh, merged, selected)
        if species_database_root:
            from bactipipe.batch_species import (
                checkm_contamination_from_summary,
                expected_organisms_from_summary,
                reconcile_batch_identifications,
                update_final_summary_verdicts,
                update_summary_identifications,
            )

            summary_path = merged / f"{run_name}_qc_summary.tsv"
            temporary_summary_path = merged / "qc_out" / "temp_qc_summary.tsv"
            expected_summary_path = (
                temporary_summary_path
                if temporary_summary_path.is_file()
                else summary_path
            )
            expected_organisms = expected_organisms_from_summary(
                expected_summary_path
            )
            reconciled = reconcile_batch_identifications(
                species_output_root=merged / "species_identification",
                database_root=species_database_root,
                expected_organisms=expected_organisms,
                platform=species_platform,
                assemblies={
                    sample_id: merged
                    / "assemblies"
                    / "genomes"
                    / f"{sample_id}.fasta"
                    for sample_id in expected_organisms
                },
                gtdb_metadata_tsv=species_gtdb_metadata,
                gtdbtk_command_prefix=(
                    "conda",
                    "run",
                    "-n",
                    species_gtdbtk_conda_env,
                    "gtdbtk",
                ),
                placement_threads=max(1, int(species_placement_threads)),
                checkm_contamination_by_sample=checkm_contamination_from_summary(
                    summary_path
                ),
                maximum_checkm_contamination=_argument_float(
                    argv, "--max_contamination", 10.0
                ),
            )
            if temporary_summary_path.is_file():
                update_summary_identifications(temporary_summary_path, reconciled)
            update_summary_identifications(summary_path, reconciled)
            if temporary_summary_path.is_file():
                update_final_summary_verdicts(
                    summary_path,
                    temporary_summary_path,
                    minimum_completeness=_argument_float(
                        argv, "--min_completeness", 90.0
                    ),
                    maximum_contamination=_argument_float(
                        argv, "--max_contamination", 10.0
                    ),
                )
        if "--defer-qc-pdf" in argv:
            (merged / f"{run_name}_qc_summary.pdf").unlink(missing_ok=True)
        _record_update(merged, fresh, run_name, selected)
        _publish_directory(merged, target)
    finally:
        shutil.rmtree(stage_parent, ignore_errors=True)


def _validate_existing_results(target: Path, run_name: str) -> None:
    if not target.is_dir() or target.is_symlink():
        raise FileNotFoundError(
            f"--update-existing requires an existing BactiPipe output directory: {target}"
        )
    summary = target / f"{run_name}_qc_summary.tsv"
    if not summary.is_file():
        raise FileNotFoundError(f"Required existing BactiPipe QC summary is missing: {summary}")


def _replace_outdir_argument(argv: list[str], fresh_parent: Path) -> list[str]:
    rewritten: list[str] = []
    skip_value = False
    for value in argv:
        if skip_value:
            skip_value = False
            continue
        if value == "--update-existing":
            continue
        if value in {"-o", "--outdir"}:
            skip_value = True
            continue
        if value.startswith("--outdir="):
            continue
        rewritten.append(value)
    rewritten.extend(("--outdir", str(fresh_parent)))
    return rewritten


def _replace_sample_outputs(fresh: Path, merged: Path, sample_name: str) -> None:
    relatives = (
        Path("qc_out") / sample_name,
        Path("qc_out") / f"{sample_name}_metrics.txt",
        Path("qc_out") / f"{sample_name}_metrics.png",
        Path("qc_out") / f"{sample_name}_quality_distribution.png",
        Path("assemblies") / sample_name,
        Path("assemblies") / "genomes" / f"{sample_name}.fasta",
        Path("assemblies") / "genomes" / f"{sample_name}.input.json",
        Path("checkM") / "bins" / f"{sample_name}.fasta",
        Path("species_identification") / sample_name,
    )
    for relative in relatives:
        source = fresh / relative
        destination = merged / relative
        if destination.is_dir():
            shutil.rmtree(destination)
        elif destination.exists():
            destination.unlink()
        if source.is_dir():
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(source, destination)
        elif source.is_file():
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)


def _merge_qc_summary(
    existing: Path,
    fresh: Path,
    merged: Path,
    run_name: str,
    sample_names: tuple[str, ...],
) -> None:
    summary_name = f"{run_name}_qc_summary.tsv"
    old_prefix, old_rows = _split_summary(
        (existing / summary_name).read_text(encoding="utf-8", errors="replace")
    )
    new_prefix, new_rows = _split_summary(
        (fresh / summary_name).read_text(encoding="utf-8", errors="replace")
    )
    selected = set(sample_names)
    replacements = {row[0]: row for row in new_rows if row and row[0] in selected}
    if set(replacements) != selected:
        missing = sorted(selected - set(replacements))
        raise ValueError("Reanalysis summary is missing sample(s): " + ", ".join(missing))
    output_rows = [replacements.get(row[0], row) if row else row for row in old_rows]
    existing_names = {row[0] for row in old_rows if row}
    output_rows.extend(replacements[name] for name in sample_names if name not in existing_names)
    with (merged / summary_name).open("w", encoding="utf-8", newline="") as handle:
        handle.write("\n".join(new_prefix or old_prefix) + "\n")
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(output_rows)


def _merge_temporary_summary(
    existing: Path,
    fresh: Path,
    merged: Path,
    sample_names: tuple[str, ...],
) -> None:
    relative = Path("qc_out") / "temp_qc_summary.tsv"
    old_path = existing / relative
    new_path = fresh / relative
    if not old_path.is_file() or not new_path.is_file():
        return
    old_prefix, old_rows = _split_temporary_summary(old_path.read_text())
    new_prefix, new_rows = _split_temporary_summary(new_path.read_text())
    selected = set(sample_names)
    replacements = {row[0]: row for row in new_rows if row and row[0] in selected}
    if set(replacements) != selected:
        missing = sorted(selected - set(replacements))
        raise ValueError(
            "Reanalysis temporary summary is missing sample(s): " + ", ".join(missing)
        )
    output_rows = [replacements.get(row[0], row) if row else row for row in old_rows]
    existing_names = {row[0] for row in old_rows if row}
    output_rows.extend(replacements[name] for name in sample_names if name not in existing_names)
    destination = merged / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\n".join(new_prefix or old_prefix) + "\n")
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(output_rows)


def _merge_failure_audit(
    existing: Path,
    fresh: Path,
    merged: Path,
    sample_names: tuple[str, ...],
) -> None:
    selected = set(sample_names)
    old_rows = _read_dict_tsv(existing / "sample_failures.tsv")
    new_rows = _read_dict_tsv(fresh / "sample_failures.tsv")
    rows = [row for row in old_rows if row.get("sample_id") not in selected]
    rows.extend(row for row in new_rows if row.get("sample_id") in selected)
    path = merged / "sample_failures.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("sample_id", "stage", "reason"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _selected_failures(path: Path, sample_names: tuple[str, ...]) -> list[dict[str, str]]:
    selected = set(sample_names)
    return [row for row in _read_dict_tsv(path) if row.get("sample_id") in selected]


def _record_update(
    root: Path,
    fresh_root: Path,
    run_name: str,
    sample_names: tuple[str, ...],
) -> None:
    completed_at = datetime.now(timezone.utc)
    audit_dir = root / "sample_updates" / completed_at.strftime("%Y%m%dT%H%M%S%fZ")
    audit_dir.mkdir(parents=True, exist_ok=True)
    for name in (f"{run_name}_qc_summary.tsv", "sample_failures.tsv"):
        source = root / name
        if source.is_file():
            shutil.copy2(source, audit_dir / name)
    selected_audit = audit_dir / "selected_reanalysis"
    selected_audit.mkdir()
    for name in (f"{run_name}_qc_summary.tsv", "sample_failures.tsv"):
        source = fresh_root / name
        if source.is_file():
            shutil.copy2(source, selected_audit / name)
    for relative in (
        Path("qc_out") / "temp_qc_summary.tsv",
        Path("species_identification") / "cross_sample_signal_review.tsv",
        Path("species_identification") / "cross_sample_signal_review.json",
    ):
        source = root / relative
        if source.is_file():
            destination = audit_dir / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    path = root / "sample_update_history.jsonl"
    with path.open("a", encoding="utf-8") as handle:
        handle.write(
            json.dumps(
                {
                    "completed_at": completed_at.isoformat(),
                    "provenance": str(audit_dir.relative_to(root)),
                    "samples": list(sample_names),
                    "status": "published",
                },
                sort_keys=True,
            )
            + "\n"
        )


def _publish_directory(staged: Path, target: Path) -> None:
    backup = target.parent / f".{target.name}.sample-update-backup-{uuid.uuid4().hex}"
    os.replace(target, backup)
    try:
        os.replace(staged, target)
    except Exception:
        os.replace(backup, target)
        raise
    shutil.rmtree(backup, ignore_errors=True)


def _split_summary(text: str) -> tuple[list[str], list[list[str]]]:
    lines = text.splitlines()
    header_index = next(
        (index for index, line in enumerate(lines) if line.startswith("Sample ID\t")),
        None,
    )
    if header_index is None:
        raise ValueError("BactiPipe QC summary header was not found.")
    rows = list(csv.reader(lines[header_index + 1 :], delimiter="\t"))
    return lines[: header_index + 1], [row for row in rows if row]


def _split_temporary_summary(text: str) -> tuple[list[str], list[list[str]]]:
    lines = text.splitlines()
    header_index = next(
        (index for index, line in enumerate(lines) if line.startswith("Sample\t")),
        None,
    )
    if header_index is None:
        raise ValueError("BactiPipe temporary QC summary header was not found.")
    rows = list(csv.reader(lines[header_index + 1 :], delimiter="\t"))
    return lines[: header_index + 1], [row for row in rows if row]


def _argument_float(argv: list[str], name: str, default: float) -> float:
    for index, value in enumerate(argv):
        if value == name and index + 1 < len(argv):
            try:
                return float(argv[index + 1])
            except ValueError:
                return default
        if value.startswith(name + "="):
            try:
                return float(value.split("=", 1)[1])
            except ValueError:
                return default
    return default


def _read_dict_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
