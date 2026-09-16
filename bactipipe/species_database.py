"""Build immutable, inactive bacterial species-identification databases."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from bactipipe.species_identification import ReferenceCatalog, ReferenceRecord


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


def build_identification_database(
    *,
    reviewed_manifest: str | Path,
    output_root: str | Path,
    threads: int = 1,
    version: str | None = None,
    skani_executable: str = "skani",
    ncbi_taxonomy_dir: str | Path,
    kraken2_build_executable: str = "kraken2-build",
    kraken_max_db_size: int | None = None,
    kraken_library_batch_size: int = 100,
) -> Path:
    """Publish immutable skani and Kraken2 indexes without activating them."""

    source_manifest = Path(reviewed_manifest).resolve()
    catalog = ReferenceCatalog.from_tsv(source_manifest)
    catalog.validate_files(require_checksums=True)
    taxonomy_source = Path(ncbi_taxonomy_dir).resolve()
    _validate_taxonomy_directory(taxonomy_source)
    if kraken_max_db_size is not None and kraken_max_db_size < 1:
        raise ValueError("kraken_max_db_size must be a positive byte count")
    if kraken_library_batch_size < 1:
        raise ValueError("kraken_library_batch_size must be at least 1")
    manifest_digest = _sha256(source_manifest)
    version = version or (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        + "-"
        + manifest_digest[:12]
    )
    if not version.replace("-", "").replace("_", "").isalnum():
        raise ValueError("Database version may contain only letters, numbers, '-' and '_'")

    versions_root = Path(output_root).resolve() / "versions"
    final = versions_root / version
    staging = versions_root / f".{version}.building"
    if final.exists() or staging.exists():
        raise FileExistsError(f"Database version already exists or is being built: {version}")
    staging.mkdir(parents=True)
    try:
        reference_root = staging / "references"
        reference_root.mkdir()
        staged_records: list[ReferenceRecord] = []
        for record in catalog.records:
            source = Path(record.fasta_path)
            source = source if source.is_absolute() else catalog.root / source
            suffix = _fasta_suffix(source)
            destination = reference_root / f"{record.accession}{suffix}"
            _copy_or_link(source, destination)
            staged_records.append(
                ReferenceRecord(
                    **{
                        **record.__dict__,
                        "fasta_path": destination.relative_to(staging).as_posix(),
                    }
                )
            )

        candidate_manifest = staging / "reference_manifest.tsv"
        _write_manifest(candidate_manifest, staged_records)
        ReferenceCatalog.from_tsv(candidate_manifest).validate_files(require_checksums=True)

        callable_paths = [
            record.fasta_path
            for record in staged_records
            if record.callable_for_identification
        ]
        if not callable_paths:
            raise ValueError("Reviewed manifest has no callable isolate references")
        list_path = staging / "skani_input_files.txt"
        list_path.write_text("".join(f"{path}\n" for path in callable_paths))
        skani_database = staging / "skani"
        command = [
            skani_executable,
            "sketch",
            "-l",
            list_path.name,
            "-o",
            skani_database.name,
            "-t",
            str(max(1, int(threads))),
        ]
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            cwd=staging,
        )
        (staging / "skani_build.stdout.log").write_text(completed.stdout or "")
        (staging / "skani_build.stderr.log").write_text(completed.stderr or "")
        if completed.returncode != 0:
            raise RuntimeError(
                f"skani sketch failed with exit {completed.returncode}; "
                f"see {staging / 'skani_build.stderr.log'}"
            )
        if not skani_database.exists():
            raise RuntimeError("skani sketch completed without producing its database")
        kraken_records = [
            record for record in staged_records if record.include_in_kraken
        ]
        if not kraken_records:
            raise ValueError("Reviewed manifest has no references selected for Kraken2")
        taxonomy_audit = staging / "ncbi_taxonomy"
        taxonomy_audit.mkdir()
        taxonomy_checksums = _copy_taxonomy(
            taxonomy_source=taxonomy_source,
            destination=taxonomy_audit,
        )
        kraken_database = staging / "kraken2"
        kraken_taxonomy = kraken_database / "taxonomy"
        kraken_taxonomy.mkdir(parents=True)
        _copy_taxonomy(
            taxonomy_source=taxonomy_audit,
            destination=kraken_taxonomy,
        )
        _add_kraken_reference_batches(
            records=kraken_records,
            database_root=staging,
            kraken_database=kraken_database,
            executable=kraken2_build_executable,
            batch_size=kraken_library_batch_size,
            stdout_path=staging / "kraken2_add.stdout.log",
            stderr_path=staging / "kraken2_add.stderr.log",
        )
        kraken_build_command = [
            kraken2_build_executable,
            "--build",
            "--threads",
            str(max(1, int(threads))),
            "--db",
            str(kraken_database),
        ]
        if kraken_max_db_size:
            kraken_build_command.extend(["--max-db-size", str(kraken_max_db_size)])
        _run_logged(
            kraken_build_command,
            stdout_path=staging / "kraken2_build.stdout.log",
            stderr_path=staging / "kraken2_build.stderr.log",
            failure_label="Kraken2 database build",
        )
        _validate_kraken_database(kraken_database)
        _run_logged(
            [kraken2_build_executable, "--clean", "--db", str(kraken_database)],
            stdout_path=staging / "kraken2_clean.stdout.log",
            stderr_path=staging / "kraken2_clean.stderr.log",
            failure_label="Kraken2 database cleanup",
        )
        _validate_kraken_database(kraken_database)

        identity = {
            "schema_version": 1,
            "name": "BactiPipe bacterial species-identification database",
            "version": version,
            "status": "candidate",
            "active": False,
            "built_at": datetime.now(timezone.utc).isoformat(),
            "reference_count": len(staged_records),
            "callable_reference_count": len(callable_paths),
            "species_count": len(
                {record.species_taxid for record in staged_records if record.callable_for_identification}
            ),
            "reporting_group_count": len(
                {
                    record.reporting_group
                    for record in staged_records
                    if record.reporting_group
                }
            ),
            "reporting_group_species_count": len(
                {
                    record.species_taxid
                    for record in staged_records
                    if record.reporting_group
                }
            ),
            "source_manifest_sha256": manifest_digest,
            "reference_manifest_sha256": _sha256(candidate_manifest),
            "skani_version": _tool_version(skani_executable),
            "kraken2_version": _tool_version(kraken2_build_executable),
            "kraken_reference_count": len(kraken_records),
            "kraken_species_count": len(
                {record.species_taxid for record in kraken_records}
            ),
            "kraken_reference_manifest_sha256": _sha256(candidate_manifest),
            "kraken_max_db_size": kraken_max_db_size,
            "kraken_library_batch_size": kraken_library_batch_size,
            "kraken_taxonomy_files": taxonomy_checksums,
            "index_relationship": (
                "The skani and Kraken2 indexes were built from the same reviewed "
                "reference manifest; Kraken2 is candidate-nomination only."
            ),
            "selection_policy": (
                "Reviewed isolate-derived type/reference/representative genomes; "
                "MAG/SAG records cannot be callable."
            ),
        }
        (staging / "database_identity.json").write_text(
            json.dumps(identity, indent=2, sort_keys=True) + "\n"
        )
        _write_checksums(staging)
        (staging / "BUILD_SUCCESS").write_text(
            f"{identity['built_at']}\t{version}\tcandidate_inactive\n"
        )
        os.replace(staging, final)
    except Exception:
        (staging / "BUILD_FAILED").write_text(
            datetime.now(timezone.utc).isoformat() + "\n"
        )
        raise
    return final


def validate_database(path: str | Path, *, allow_candidate: bool = False) -> dict[str, object]:
    root = Path(path)
    identity_path = root / "database_identity.json"
    success_path = root / "BUILD_SUCCESS"
    manifest_path = root / "reference_manifest.tsv"
    skani_path = root / "skani"
    kraken_path = root / "kraken2"
    for required in (
        identity_path,
        success_path,
        manifest_path,
        skani_path,
        kraken_path,
    ):
        if not required.exists():
            raise ValueError(f"Species-identification database is incomplete: {required}")
    identity = json.loads(identity_path.read_text())
    if not allow_candidate and not (
        identity.get("status") == "production" and identity.get("active") is True
    ):
        raise ValueError(
            f"Database {identity.get('version', 'unknown')} is not an active production database"
        )
    catalog = ReferenceCatalog.from_tsv(manifest_path)
    catalog.validate_files(require_checksums=True)
    _validate_kraken_database(kraken_path)
    if identity.get("reference_manifest_sha256") != _sha256(manifest_path):
        raise ValueError("Reference manifest checksum does not match database identity")
    if identity.get("kraken_reference_manifest_sha256") != _sha256(manifest_path):
        raise ValueError("Kraken2 index manifest checksum does not match database identity")
    return identity


def _write_manifest(path: Path, records: list[ReferenceRecord]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=MANIFEST_COLUMNS)
        writer.writeheader()
        for record in records:
            row = record.__dict__.copy()
            for field in (
                "type_material",
                "canonical_for_tree",
                "callable_for_identification",
                "include_in_kraken",
            ):
                row[field] = "yes" if row[field] else "no"
            writer.writerow(row)


def _copy_or_link(source: Path, destination: Path) -> None:
    try:
        os.link(source, destination)
    except OSError:
        shutil.copy2(source, destination)


def _fasta_suffix(path: Path) -> str:
    """Return only the biological file suffix, excluding accession versions."""

    lower_name = path.name.casefold()
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if lower_name.endswith(suffix):
            return path.name[-len(suffix) :]
    return ".fna"


def _validate_taxonomy_directory(path: Path) -> None:
    for name in ("nodes.dmp", "names.dmp"):
        candidate = path / name
        if not candidate.is_file() or candidate.stat().st_size == 0:
            raise ValueError(f"NCBI taxonomy file is missing or empty: {candidate}")


def _copy_taxonomy(*, taxonomy_source: Path, destination: Path) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for name in ("nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"):
        source = taxonomy_source / name
        if not source.is_file():
            continue
        target = destination / name
        shutil.copy2(source, target)
        checksums[name] = _sha256(target)
    return checksums


def _write_kraken_library(
    destination: Path,
    *,
    database_root: Path,
    records: list[ReferenceRecord],
) -> None:
    """Write one explicitly taxid-labelled library from the reviewed manifest."""

    with destination.open("w") as output:
        for record in records:
            source = database_root / record.fasta_path
            opener = gzip.open if source.suffix == ".gz" else open
            sequence_number = 0
            with opener(source, "rt") as handle:
                for line in handle:
                    if line.startswith(">"):
                        sequence_number += 1
                        original = line[1:].strip()
                        sequence_id = (
                            f"{record.accession}_{sequence_number}|kraken:taxid|"
                            f"{record.species_taxid}"
                        )
                        output.write(f">{sequence_id} {original}\n")
                    else:
                        output.write(line)
            if sequence_number == 0:
                raise ValueError(f"Reference FASTA contains no records: {source}")


def _add_kraken_reference_batches(
    *,
    records: list[ReferenceRecord],
    database_root: Path,
    kraken_database: Path,
    executable: str,
    batch_size: int,
    stdout_path: Path,
    stderr_path: Path,
) -> None:
    """Import bounded temporary libraries to avoid a full extra FASTA copy."""

    batch_root = database_root / ".kraken-input-batches"
    batch_root.mkdir()
    stdout_path.write_text("")
    stderr_path.write_text("")
    for offset in range(0, len(records), batch_size):
        batch_number = offset // batch_size + 1
        batch = records[offset : offset + batch_size]
        batch_path = batch_root / f"batch_{batch_number:06d}.fna"
        _write_kraken_library(
            batch_path,
            database_root=database_root,
            records=batch,
        )
        _run_logged(
            [
                executable,
                "--add-to-library",
                str(batch_path),
                "--db",
                str(kraken_database),
            ],
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            failure_label=f"Kraken2 library import batch {batch_number}",
            append=True,
        )
        batch_path.unlink()
    batch_root.rmdir()


def _run_logged(
    command: list[str],
    *,
    stdout_path: Path,
    stderr_path: Path,
    failure_label: str,
    append: bool = False,
) -> None:
    completed = subprocess.run(command, text=True, capture_output=True)
    mode = "a" if append else "w"
    with stdout_path.open(mode) as stdout:
        stdout.write(completed.stdout or "")
    with stderr_path.open(mode) as stderr:
        stderr.write(completed.stderr or "")
    if completed.returncode != 0:
        raise RuntimeError(
            f"{failure_label} failed with exit {completed.returncode}; "
            f"see {stderr_path}"
        )


def _validate_kraken_database(path: Path) -> None:
    for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        candidate = path / name
        if not candidate.is_file() or candidate.stat().st_size == 0:
            raise ValueError(f"Kraken2 database is incomplete: {candidate}")


def _tool_version(executable: str) -> str:
    completed = subprocess.run(
        [executable, "--version"], text=True, capture_output=True, check=False
    )
    if completed.returncode != 0:
        return "unknown"
    return (completed.stdout or completed.stderr or "unknown").strip().splitlines()[0]


def _write_checksums(root: Path) -> None:
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name not in {"checksums.sha256", "BUILD_SUCCESS"}:
            rows.append(f"{_sha256(path)}  {path.relative_to(root).as_posix()}\n")
    (root / "checksums.sha256").write_text("".join(rows))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
