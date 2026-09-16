"""Checkpointed retrieval of current NCBI assembly reports."""

from __future__ import annotations

import hashlib
import csv
import json
import os
import re
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


ASSEMBLY_ACCESSION = re.compile(r"GC[AF]_\d+\.\d+")


def read_assembly_accessions(path: str | Path) -> tuple[str, ...]:
    accessions: list[str] = []
    seen: set[str] = set()
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            accession = line.strip()
            if not accession or accession.startswith("#"):
                continue
            if not ASSEMBLY_ACCESSION.fullmatch(accession):
                raise ValueError(f"Invalid assembly accession on line {line_number}: {accession}")
            if accession not in seen:
                accessions.append(accession)
                seen.add(accession)
    if not accessions:
        raise ValueError("Assembly accession input is empty")
    return tuple(accessions)


def fetch_assembly_reports(
    accessions_path: str | Path,
    output_root: str | Path,
    *,
    datasets_executable: str = "datasets",
    batch_size: int = 250,
    retries: int = 4,
) -> dict[str, object]:
    """Fetch reports in independently resumable batches and publish an audit."""

    if batch_size < 1:
        raise ValueError("batch_size must be at least 1")
    if retries < 1:
        raise ValueError("retries must be at least 1")
    accessions_source = Path(accessions_path).resolve()
    requested = read_assembly_accessions(accessions_source)
    output = Path(output_root).resolve()
    batches_root = output / "batches"
    batches_root.mkdir(parents=True, exist_ok=True)
    all_reports: dict[str, dict[str, object]] = {}
    accession_resolution: dict[str, str] = {}
    returned_accessions: set[str] = set()

    for index, offset in enumerate(range(0, len(requested), batch_size), start=1):
        batch = requested[offset : offset + batch_size]
        batch_root = batches_root / f"{index:05d}"
        batch_root.mkdir(exist_ok=True)
        accession_file = batch_root / "accessions.txt"
        report_file = batch_root / "assembly_data_report.jsonl"
        success_file = batch_root / "BATCH_SUCCESS.json"
        expected_accession_text = "".join(f"{accession}\n" for accession in batch)
        if accession_file.exists() and accession_file.read_text() != expected_accession_text:
            raise ValueError(f"Batch checkpoint input changed: {batch_root}")
        _atomic_text(accession_file, expected_accession_text)

        batch_reports: dict[str, dict[str, object]]
        checkpointed = success_file.exists() and report_file.exists()
        if checkpointed:
            checkpoint = json.loads(success_file.read_text())
            if checkpoint.get("requested_sha256") != _sha256(accession_file):
                raise ValueError(f"Batch checkpoint checksum mismatch: {batch_root}")
            if checkpoint.get("report_sha256") != _sha256(report_file):
                raise ValueError(f"Batch report checksum mismatch: {batch_root}")
            batch_reports = _read_reports(report_file)
        else:
            batch_reports = _fetch_batch(
                accession_file,
                report_file,
                datasets_executable=datasets_executable,
                retries=retries,
                stderr_path=batch_root / "datasets.stderr.log",
            )
        resolution = _resolve_batch_reports(batch, batch_reports, batch_root=batch_root)
        missing = [accession for accession in batch if accession not in resolution]
        if not checkpointed:
            _atomic_text(
                success_file,
                json.dumps(
                    {
                        "requested": len(batch),
                        "retrieved": len(resolution),
                        "missing": missing,
                        "accession_replacements": {
                            requested_accession: returned_accession
                            for requested_accession, returned_accession in resolution.items()
                            if requested_accession != returned_accession
                        },
                        "requested_sha256": _sha256(accession_file),
                        "report_sha256": _sha256(report_file),
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
            )
        overlap = returned_accessions.intersection(resolution.values())
        if overlap:
            raise ValueError(f"NCBI returned duplicate accessions across batches: {sorted(overlap)[:3]}")
        returned_accessions.update(resolution.values())
        for requested_accession, returned_accession in resolution.items():
            all_reports[requested_accession] = batch_reports[returned_accession]
            accession_resolution[requested_accession] = returned_accession

    retrieved = [accession for accession in requested if accession in all_reports]
    missing = [accession for accession in requested if accession not in all_reports]
    combined = output / "assembly_data_report.jsonl"
    _atomic_text(
        combined,
        "".join(
            json.dumps(all_reports[accession], sort_keys=True, separators=(",", ":")) + "\n"
            for accession in retrieved
        ),
    )
    _atomic_text(output / "retrieved_accessions.txt", "".join(f"{item}\n" for item in retrieved))
    _atomic_text(output / "missing_accessions.txt", "".join(f"{item}\n" for item in missing))
    _atomic_tsv(
        output / "accession_resolution.tsv",
        ("requested_accession", "returned_accession", "status"),
        (
            {
                "requested_accession": accession,
                "returned_accession": accession_resolution.get(accession, ""),
                "status": (
                    "missing"
                    if accession not in accession_resolution
                    else "current_version_substituted"
                    if accession_resolution[accession] != accession
                    else "exact"
                ),
            }
            for accession in requested
        ),
    )
    identity: dict[str, object] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "datasets_version": _datasets_version(datasets_executable),
        "batch_size": batch_size,
        "requested_count": len(requested),
        "retrieved_count": len(retrieved),
        "missing_count": len(missing),
        "current_version_substitution_count": sum(
            requested_accession != returned_accession
            for requested_accession, returned_accession in accession_resolution.items()
        ),
        "accession_input_sha256": _sha256(accessions_source),
        "assembly_data_report_sha256": _sha256(combined),
        "status": "success" if not missing else "success_with_missing_accessions",
    }
    _atomic_text(
        output / "fetch_summary.json",
        json.dumps(identity, indent=2, sort_keys=True) + "\n",
    )
    _atomic_text(
        output / "FETCH_SUCCESS",
        f"{identity['created_at']}\t{identity['status']}\n",
    )
    return identity


def _fetch_batch(
    accessions_path: Path,
    report_path: Path,
    *,
    datasets_executable: str,
    retries: int,
    stderr_path: Path,
) -> dict[str, dict[str, object]]:
    temporary = report_path.with_suffix(report_path.suffix + ".tmp")
    command = [
        datasets_executable,
        "summary",
        "genome",
        "accession",
        "--inputfile",
        str(accessions_path),
        "--as-json-lines",
        "--assembly-version",
        "current",
        "--mag",
        "exclude",
        "--exclude-atypical",
        "--exclude-multi-isolate",
    ]
    last_error = ""
    for attempt in range(1, retries + 1):
        with temporary.open("w") as stdout, stderr_path.open("a") as stderr:
            completed = subprocess.run(command, text=True, stdout=stdout, stderr=stderr)
        if completed.returncode == 0:
            reports = _read_reports(temporary)
            os.replace(temporary, report_path)
            return reports
        last_error = f"exit {completed.returncode}"
        temporary.unlink(missing_ok=True)
        if attempt < retries:
            time.sleep(min(30, 2**attempt))
    raise RuntimeError(
        f"NCBI Datasets report retrieval failed after {retries} attempts ({last_error}); "
        f"see {stderr_path}"
    )


def _read_reports(path: Path) -> dict[str, dict[str, object]]:
    reports: dict[str, dict[str, object]] = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            value = json.loads(line)
            if not isinstance(value, dict):
                raise ValueError(f"NCBI report line {line_number} is not an object")
            accession = str(value.get("accession") or "").strip()
            if not ASSEMBLY_ACCESSION.fullmatch(accession):
                raise ValueError(f"NCBI report line {line_number} has invalid accession")
            if accession in reports:
                raise ValueError(f"NCBI report contains duplicate accession {accession}")
            reports[accession] = value
    return reports


def _resolve_batch_reports(
    requested: Iterable[str],
    reports: dict[str, dict[str, object]],
    *,
    batch_root: Path,
) -> dict[str, str]:
    requested = tuple(requested)
    by_stem: dict[str, list[str]] = {}
    for returned in reports:
        by_stem.setdefault(returned.rsplit(".", 1)[0], []).append(returned)
    resolution: dict[str, str] = {}
    requested_stems = {accession.rsplit(".", 1)[0] for accession in requested}
    for accession in requested:
        matches = by_stem.get(accession.rsplit(".", 1)[0], [])
        current_matches = [
            returned
            for returned in matches
            if _dict(reports[returned].get("assembly_info")).get("assembly_status")
            in {None, "", "current"}
        ]
        if len(current_matches) == 1:
            resolution[accession] = current_matches[0]
        elif len(current_matches) > 1:
            raise ValueError(
                f"NCBI returned multiple current versions for {accession}: {current_matches}"
            )
        elif len(matches) == 1:
            # Retain the record so reconciliation can explicitly reject a
            # suppressed assembly instead of silently treating it as absent.
            resolution[accession] = matches[0]
    unexpected = set(reports).difference(resolution.values())
    truly_unexpected = {
        accession
        for accession in unexpected
        if accession.rsplit(".", 1)[0] not in requested_stems
    }
    if truly_unexpected:
        raise ValueError(
            f"NCBI returned unrequested accessions for {batch_root}: "
            f"{sorted(truly_unexpected)[:3]}"
        )
    return resolution


def _dict(value: object) -> dict[str, object]:
    return value if isinstance(value, dict) else {}


def _datasets_version(executable: str) -> str:
    completed = subprocess.run(
        [executable, "version"], text=True, capture_output=True, check=False
    )
    text = " ".join((completed.stdout or completed.stderr or "unknown").split())
    return text[:200]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_text(path: str | Path, content: str) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.write_text(content)
    os.replace(temporary, destination)


def _atomic_tsv(
    path: str | Path,
    fields: tuple[str, ...],
    rows: Iterable[dict[str, object]],
) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, destination)
