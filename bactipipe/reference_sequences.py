"""Checkpointed NCBI genome download and reviewed-manifest generation."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import subprocess
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Sequence

from bactipipe.reporting_groups import ReportingGroupCatalog


ASSEMBLY_ACCESSION = re.compile(r"GC[AF]_\d+\.\d+")
MANIFEST_COLUMNS = [
    "accession",
    "species_taxid",
    "species_name",
    "fasta_path",
    "source",
    "reference_role",
    "assembly_level",
    "sequence_origin",
    "type_material",
    "canonical_for_tree",
    "callable_for_identification",
    "include_in_kraken",
    "reporting_group",
    "ncbi_taxonomy_name",
    "gtdb_taxonomy",
    "sha256",
]


class GenomeDownloadError(RuntimeError):
    pass


def download_reference_genomes(
    inventory_tsv: str | Path,
    output_root: str | Path,
    *,
    datasets_executable: str = "datasets",
    batch_size: int = 50,
    retries: int = 4,
    retry_delay_seconds: float = 2.0,
    timeout_seconds: float = 1800.0,
    runner: Callable[..., subprocess.CompletedProcess[str]] | None = None,
    sleeper: Callable[[float], None] = time.sleep,
) -> dict[str, object]:
    """Download one independently verified FASTA per assembly accession."""

    if batch_size < 1 or retries < 1:
        raise ValueError("batch_size and retries must be at least 1")
    inventory_path = Path(inventory_tsv).resolve()
    rows = _read_inventory(inventory_path)
    accessions = tuple(row["accession"] for row in rows)
    inventory_sha256 = _sha256(inventory_path)
    output = Path(output_root).resolve()
    references = output / "references"
    batches = output / "batches"
    references.mkdir(parents=True, exist_ok=True)
    batches.mkdir(parents=True, exist_ok=True)
    state_path = output / "download_state.json"
    state = _load_state(
        state_path,
        inventory_path=inventory_path,
        inventory_sha256=inventory_sha256,
        accession_count=len(accessions),
    )
    state.update(status="running", updated_at=_utc_now())
    _atomic_json(state_path, state)

    run = runner or _run
    results = []
    try:
        for number, offset in enumerate(range(0, len(accessions), batch_size), start=1):
            results.append(
                _download_batch(
                    f"batch-{number:06d}",
                    accessions[offset : offset + batch_size],
                    references=references,
                    batches=batches,
                    datasets_executable=datasets_executable,
                    retries=retries,
                    retry_delay_seconds=retry_delay_seconds,
                    timeout_seconds=timeout_seconds,
                    runner=run,
                    sleeper=sleeper,
                )
            )
    except BaseException as exc:
        state.update(
            status="interrupted" if isinstance(exc, KeyboardInterrupt) else "incomplete",
            last_error=str(exc),
            updated_at=_utc_now(),
        )
        _atomic_json(state_path, state)
        raise

    downloaded = {
        accession: record
        for result in results
        for accession, record in result["downloaded"].items()
    }
    missing = {
        accession for result in results for accession in result["missing_accessions"]
    }
    if set(downloaded).intersection(missing):
        raise GenomeDownloadError("An accession was recorded as both downloaded and missing")
    unaccounted = set(accessions).difference(downloaded, missing)
    if unaccounted:
        raise GenomeDownloadError(
            f"Download left {len(unaccounted)} accessions without an outcome"
        )

    audit_rows = []
    for accession in accessions:
        record = downloaded.get(accession)
        audit_rows.append(
            {
                "accession": accession,
                "status": "downloaded" if record else "missing",
                "fasta_path": record["fasta_path"] if record else "",
                "sha256": record["sha256"] if record else "",
                "sequence_count": record["sequence_count"] if record else "",
                "total_bases": record["total_bases"] if record else "",
                "batch_id": record["batch_id"] if record else "",
            }
        )
    audit_path = output / "sequence_download_audit.tsv"
    _atomic_tsv(
        audit_path,
        [
            "accession",
            "status",
            "fasta_path",
            "sha256",
            "sequence_count",
            "total_bases",
            "batch_id",
        ],
        audit_rows,
    )
    _atomic_text(
        output / "downloaded_accessions.txt",
        "".join(f"{accession}\n" for accession in accessions if accession in downloaded),
    )
    _atomic_text(
        output / "missing_accessions.txt",
        "".join(f"{accession}\n" for accession in accessions if accession in missing),
    )
    summary: dict[str, object] = {
        "schema_version": 1,
        "status": "success" if not missing else "success_with_missing_accessions",
        "completed_at": _utc_now(),
        "inventory": str(inventory_path),
        "inventory_sha256": inventory_sha256,
        "datasets_version": _tool_version(datasets_executable),
        "batch_size": batch_size,
        "requested_count": len(accessions),
        "downloaded_count": len(downloaded),
        "missing_count": len(missing),
        "audit_sha256": _sha256(audit_path),
    }
    _atomic_json(output / "download_summary.json", summary)
    state.update(
        status=summary["status"],
        completed_at=summary["completed_at"],
        updated_at=summary["completed_at"],
        last_error="",
    )
    _atomic_json(state_path, state)
    _atomic_text(
        output / "DOWNLOAD_SUCCESS",
        f"{summary['completed_at']}\t{summary['status']}\n",
    )
    return summary


def write_reviewed_manifest(
    inventory_tsv: str | Path,
    sequence_audit_tsv: str | Path,
    output_path: str | Path,
) -> dict[str, int]:
    """Join reviewed metadata to downloaded FASTAs without changing taxonomy."""

    inventory_path = Path(inventory_tsv).resolve()
    audit_path = Path(sequence_audit_tsv).resolve()
    destination = Path(output_path).resolve()
    inventory = _read_inventory(inventory_path)
    reporting_groups = ReportingGroupCatalog.from_package()
    audit = {
        row["accession"]: row
        for row in _read_tsv(audit_path)
        if row.get("status") == "downloaded"
    }
    manifest_rows = []
    missing = []
    canonical_species: set[str] = set()
    for row in inventory:
        accession = row["accession"]
        sequence = audit.get(accession)
        if sequence is None:
            missing.append(accession)
            continue
        expected_length_text = row.get("current_total_length") or row.get("total_length")
        if expected_length_text:
            try:
                expected_length = int(expected_length_text)
            except ValueError as exc:
                raise ValueError(
                    f"Inventory accession {accession} has an invalid expected length"
                ) from exc
            if int(sequence["total_bases"]) != expected_length:
                raise ValueError(
                    f"Downloaded FASTA length for {accession} is {sequence['total_bases']}; "
                    f"current NCBI metadata requires {expected_length}"
                )
        species_name = _manifest_species_name(row)
        canonical = row.get("reference_role") == "canonical"
        if canonical:
            if row["species_taxid"] in canonical_species:
                raise ValueError(
                    f"More than one canonical reference for taxid {row['species_taxid']}"
                )
            canonical_species.add(row["species_taxid"])
        source_fasta = audit_path.parent / sequence["fasta_path"]
        relative_fasta = os.path.relpath(source_fasta, start=destination.parent)
        manifest_rows.append(
            {
                "accession": accession,
                "species_taxid": row["species_taxid"],
                "species_name": species_name,
                "fasta_path": Path(relative_fasta).as_posix(),
                "source": row.get("current_source") or row.get("source", ""),
                "reference_role": _manifest_reference_role(row),
                "assembly_level": row.get("current_assembly_level")
                or row.get("assembly_level", ""),
                "sequence_origin": "isolate",
                "type_material": _yes_no(
                    row.get("current_type_material") or row.get("type_material")
                ),
                "canonical_for_tree": "yes" if canonical else "no",
                "callable_for_identification": "yes",
                "include_in_kraken": "yes" if canonical else "no",
                "reporting_group": reporting_groups.display_name_for_taxid(
                    row["species_taxid"]
                ),
                "ncbi_taxonomy_name": row.get("current_submitted_species")
                or row.get("ncbi_species_name", ""),
                "gtdb_taxonomy": row.get("gtdb_species_name", ""),
                "sha256": sequence["sha256"],
            }
        )
    if missing:
        raise ValueError(
            f"Cannot write a reviewed manifest; {len(missing)} assembly FASTAs are missing"
        )
    _atomic_tsv(destination, MANIFEST_COLUMNS, manifest_rows)
    summary = {
        "schema_version": 1,
        "created_at": _utc_now(),
        "inventory_sha256": _sha256(inventory_path),
        "sequence_audit_sha256": _sha256(audit_path),
        "manifest_sha256": _sha256(destination),
        "inventory_references": len(inventory),
        "manifest_references": len(manifest_rows),
        "manifest_species": len({row["species_taxid"] for row in manifest_rows}),
        "reporting_group_members": sum(
            bool(row["reporting_group"]) for row in manifest_rows
        ),
        "reporting_groups_sha256": reporting_groups.checksum,
        "missing_sequences": len(missing),
    }
    _atomic_json(destination.with_suffix(".summary.json"), summary)
    return summary


def _download_batch(
    batch_id: str,
    accessions: Sequence[str],
    *,
    references: Path,
    batches: Path,
    datasets_executable: str,
    retries: int,
    retry_delay_seconds: float,
    timeout_seconds: float,
    runner: Callable[..., subprocess.CompletedProcess[str]],
    sleeper: Callable[[float], None],
) -> dict[str, object]:
    request_digest = _line_digest(accessions)
    status_path = batches / f"{batch_id}.json"
    cached = _cached_batch(status_path, references, request_digest, len(accessions))
    if cached is not None:
        return cached
    request_path = batches / f"{batch_id}.accessions.txt"
    _atomic_text(request_path, "".join(f"{accession}\n" for accession in accessions))
    archive_path = batches / f".{batch_id}.part.zip"
    last_error = ""
    for attempt in range(1, retries + 1):
        archive_path.unlink(missing_ok=True)
        command = [
            datasets_executable,
            "download",
            "genome",
            "accession",
            "--inputfile",
            str(request_path),
            "--include",
            "genome",
            "--filename",
            str(archive_path),
            "--no-progressbar",
        ]
        try:
            completed = runner(command, timeout=timeout_seconds)
        except (OSError, subprocess.SubprocessError) as exc:
            last_error = str(exc)
            if attempt < retries:
                sleeper(min(60.0, retry_delay_seconds * (2 ** (attempt - 1))))
            continue
        if completed.returncode == 0 and archive_path.is_file():
            try:
                extracted = _extract_archive(archive_path, accessions, references, batch_id)
                archive_path.unlink(missing_ok=True)
                result = {
                    "schema_version": 1,
                    "status": "complete",
                    "batch_id": batch_id,
                    "request_sha256": request_digest,
                    "requested_count": len(accessions),
                    "downloaded": extracted,
                    "missing_accessions": [
                        accession for accession in accessions if accession not in extracted
                    ],
                    "completed_at": _utc_now(),
                }
                _atomic_json(status_path, result)
                return result
            except (GenomeDownloadError, OSError, zipfile.BadZipFile) as exc:
                last_error = str(exc)
        else:
            last_error = (
                completed.stderr or completed.stdout or f"exit {completed.returncode}"
            )[-2000:]
        if attempt < retries:
            sleeper(min(60.0, retry_delay_seconds * (2 ** (attempt - 1))))
    archive_path.unlink(missing_ok=True)
    if len(accessions) > 1:
        midpoint = len(accessions) // 2
        left = _download_batch(
            batch_id + "-a",
            accessions[:midpoint],
            references=references,
            batches=batches,
            datasets_executable=datasets_executable,
            retries=retries,
            retry_delay_seconds=retry_delay_seconds,
            timeout_seconds=timeout_seconds,
            runner=runner,
            sleeper=sleeper,
        )
        right = _download_batch(
            batch_id + "-b",
            accessions[midpoint:],
            references=references,
            batches=batches,
            datasets_executable=datasets_executable,
            retries=retries,
            retry_delay_seconds=retry_delay_seconds,
            timeout_seconds=timeout_seconds,
            runner=runner,
            sleeper=sleeper,
        )
        downloaded = {**left["downloaded"], **right["downloaded"]}
        result = {
            "schema_version": 1,
            "status": "complete",
            "batch_id": batch_id,
            "request_sha256": request_digest,
            "requested_count": len(accessions),
            "downloaded": downloaded,
            "missing_accessions": [
                accession for accession in accessions if accession not in downloaded
            ],
            "recovered_by_split": True,
            "completed_at": _utc_now(),
        }
        _atomic_json(status_path, result)
        return result
    raise GenomeDownloadError(f"Failed to download {accessions[0]}: {last_error}")


def _extract_archive(
    archive_path: Path,
    requested: Sequence[str],
    references: Path,
    batch_id: str,
) -> dict[str, dict[str, object]]:
    requested_set = set(requested)
    members: dict[str, str] = {}
    with zipfile.ZipFile(archive_path) as archive:
        for name in archive.namelist():
            parts = Path(name).parts
            if not parts or not (
                parts[-1] == "genomic.fna" or parts[-1].endswith("_genomic.fna")
            ):
                continue
            accession = next((part for part in parts if ASSEMBLY_ACCESSION.fullmatch(part)), "")
            if accession in requested_set:
                if accession in members:
                    raise GenomeDownloadError(f"Archive has duplicate FASTA for {accession}")
                members[accession] = name
        extracted: dict[str, dict[str, object]] = {}
        for accession in requested:
            member = members.get(accession)
            if not member:
                continue
            data = archive.read(member)
            sequence_count, total_bases = _validate_fasta(data, accession)
            destination = references / f"{accession}.fna"
            temporary = destination.with_suffix(".fna.tmp")
            temporary.write_bytes(data)
            os.replace(temporary, destination)
            extracted[accession] = {
                "fasta_path": destination.relative_to(references.parent).as_posix(),
                "sha256": _sha256(destination),
                "sequence_count": sequence_count,
                "total_bases": total_bases,
                "batch_id": batch_id,
            }
    return extracted


def _validate_fasta(data: bytes, accession: str) -> tuple[int, int]:
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError as exc:
        raise GenomeDownloadError(f"FASTA for {accession} is not ASCII") from exc
    count = 0
    bases = 0
    header_seen = False
    for line in text.splitlines():
        if line.startswith(">"):
            count += 1
            header_seen = True
        elif line.strip():
            if not header_seen:
                raise GenomeDownloadError(f"FASTA for {accession} starts before a header")
            bases += len("".join(line.split()))
    if count == 0 or bases == 0:
        raise GenomeDownloadError(f"FASTA for {accession} is empty")
    return count, bases


def _read_inventory(path: Path) -> list[dict[str, str]]:
    rows = _read_tsv(path)
    seen = set()
    for row in rows:
        accession = row.get("accession", "")
        if not ASSEMBLY_ACCESSION.fullmatch(accession) or accession in seen:
            raise ValueError(f"Invalid or duplicate inventory accession: {accession}")
        if not str(row.get("species_taxid", "")).isdigit():
            raise ValueError(f"Inventory accession {accession} lacks a numeric species taxid")
        seen.add(accession)
    if not rows:
        raise ValueError("Reference inventory is empty")
    return rows


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _load_state(
    path: Path,
    *,
    inventory_path: Path,
    inventory_sha256: str,
    accession_count: int,
) -> dict[str, object]:
    if path.exists():
        state = json.loads(path.read_text())
        if state.get("inventory_sha256") != inventory_sha256:
            raise ValueError("Download checkpoint inventory changed")
        return state
    return {
        "schema_version": 1,
        "inventory": str(inventory_path),
        "inventory_sha256": inventory_sha256,
        "accession_count": accession_count,
        "created_at": _utc_now(),
        "status": "pending",
    }


def _cached_batch(
    path: Path, references: Path, request_sha256: str, requested_count: int
) -> dict[str, object] | None:
    if not path.is_file():
        return None
    try:
        result = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if (
        result.get("status") != "complete"
        or result.get("request_sha256") != request_sha256
        or result.get("requested_count") != requested_count
    ):
        return None
    for record in result.get("downloaded", {}).values():
        fasta = references.parent / record["fasta_path"]
        if not fasta.is_file() or _sha256(fasta) != record["sha256"]:
            return None
    return result


def _manifest_species_name(row: dict[str, str]) -> str:
    source_taxid = row.get("source_taxid", "")
    source_species_taxid = row.get("source_species_taxid", "")
    controlled_names = [
        name.strip()
        for name in row.get("controlled_catalog_name", "").split(";")
        if name.strip()
    ]
    if row.get("controlled_match") and len(controlled_names) == 1:
        return controlled_names[0]
    if (
        row.get("controlled_match", "").endswith("ncbi_taxid")
        and source_taxid
        and source_taxid != source_species_taxid
        and controlled_names
    ):
        return controlled_names[0]
    return (
        row.get("current_submitted_species")
        or row.get("ncbi_species_name")
        or (controlled_names[0] if controlled_names else "")
    )


def _manifest_reference_role(row: dict[str, str]) -> str:
    if row.get("reference_role") != "canonical":
        return "diversity"
    if _yes_no(row.get("current_type_material") or row.get("type_material")) == "yes":
        return "type"
    category = (
        row.get("current_refseq_category") or row.get("refseq_category", "")
    ).casefold()
    if category == "reference genome":
        return "reference"
    return "representative"


def _yes_no(value: object) -> str:
    return "yes" if str(value or "").casefold() in {"1", "true", "t", "yes", "y"} else "no"


def _run(command: list[str], *, timeout: float) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, timeout=timeout)


def _tool_version(executable: str) -> str:
    result = subprocess.run(
        [executable, "version"], text=True, capture_output=True, check=False
    )
    return " ".join((result.stdout or result.stderr or "unknown").split())[:200]


def _line_digest(values: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(value.encode())
        digest.update(b"\n")
    return digest.hexdigest()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _atomic_text(path: Path, text: str) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text)
    os.replace(temporary, path)


def _atomic_json(path: Path, value: object) -> None:
    _atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def _atomic_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)
