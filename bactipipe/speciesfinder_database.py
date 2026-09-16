"""Immutable preparation and selective reference extraction for SpeciesFinder."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version as package_version
from pathlib import Path
from typing import Sequence

from bactipipe.reporting_groups import ReportingGroupCatalog
from bactipipe.species_database import MANIFEST_COLUMNS
from bactipipe.speciesfinder_evidence import SpeciesFinderHit, select_candidate_hits
from bactipipe.speciesfinder_evidence import (
    is_excluded_reference_source,
    is_reportable_species_name,
)


DATABASE_KIND = "speciesfinder_kma"
REQUIRED_KMA_SUFFIXES = (".comp.b", ".seq.b", ".length.b", ".name")


@dataclass(frozen=True)
class SpeciesFinderDatabase:
    root: Path
    database_prefix: Path
    taxonomy_file: Path
    template_lookup: Path
    identity: dict[str, object]


@dataclass(frozen=True)
class ExtractedReference:
    accession: str
    species_taxid: str
    species_name: str
    fasta_path: Path
    source_template_accession: str
    template_id: int
    reporting_group: str


def prepare_speciesfinder_database(
    *,
    database_prefix: str | Path,
    taxonomy_file: str | Path,
    output_root: str | Path,
    archive_sha256_file: str | Path,
    version: str | None = None,
    status: str = "candidate",
) -> Path:
    """Hard-link a verified bacterial KMA database into an immutable version."""

    source_prefix = Path(database_prefix).resolve()
    source_taxonomy = Path(taxonomy_file).resolve()
    archive_checksum_path = Path(archive_sha256_file).resolve()
    _validate_source(source_prefix, source_taxonomy)
    archive_sha256 = _read_sha256(archive_checksum_path)
    if status not in {"candidate", "production"}:
        raise ValueError("SpeciesFinder database status must be candidate or production")
    package = _speciesfinder_version()
    version = version or (
        f"SpeciesFinder-bacteria-v{package}-{archive_sha256[:12]}"
    )
    if not version.replace("-", "").replace("_", "").replace(".", "").isalnum():
        raise ValueError("Database version contains unsupported characters")

    versions = Path(output_root).resolve() / "versions"
    final = versions / version
    staging = versions / f".{version}.building-{uuid.uuid4().hex}"
    if final.exists():
        raise FileExistsError(f"SpeciesFinder database version already exists: {final}")
    staging.mkdir(parents=True)
    try:
        database_dir = staging / "bacteria"
        database_dir.mkdir()
        staged_prefix = database_dir / source_prefix.name
        component_metadata: dict[str, dict[str, object]] = {}
        for source in sorted(source_prefix.parent.glob(source_prefix.name + ".*")):
            if not source.is_file():
                continue
            destination = database_dir / source.name
            _copy_or_link(source, destination)
            component_metadata[destination.relative_to(staging).as_posix()] = {
                "size_bytes": destination.stat().st_size,
                "sha256": _sha256(destination),
            }
        staged_taxonomy = database_dir / source_taxonomy.name
        if not staged_taxonomy.exists():
            _copy_or_link(source_taxonomy, staged_taxonomy)
            component_metadata[staged_taxonomy.relative_to(staging).as_posix()] = {
                "size_bytes": staged_taxonomy.stat().st_size,
                "sha256": _sha256(staged_taxonomy),
            }
        _validate_source(staged_prefix, staged_taxonomy)

        taxonomy_audit = staging / "taxonomy_source_audit.tsv"
        taxonomy_counts = audit_speciesfinder_taxonomy(
            staged_taxonomy, taxonomy_audit
        )
        component_metadata[taxonomy_audit.relative_to(staging).as_posix()] = {
            "size_bytes": taxonomy_audit.stat().st_size,
            "sha256": _sha256(taxonomy_audit),
        }

        lookup = staging / "template_lookup.sqlite3"
        template_count, accession_count = _build_template_lookup(
            Path(f"{staged_prefix}.name"), lookup
        )
        component_metadata[lookup.relative_to(staging).as_posix()] = {
            "size_bytes": lookup.stat().st_size,
            "sha256": _sha256(lookup),
        }
        built_at = datetime.now(timezone.utc).isoformat()
        identity = {
            "schema_version": 1,
            "database_kind": DATABASE_KIND,
            "name": "SpeciesFinder bacterial candidate database",
            "version": version,
            "status": status,
            "active": status == "production",
            "built_at": built_at,
            "speciesfinder_version": package,
            "kma_version": _tool_version(("kma", "-v")),
            "source_url": "https://download.genepi.dk/speciesfinder_db.tar.gz",
            "source_archive_sha256": archive_sha256,
            "database_prefix": staged_prefix.relative_to(staging).as_posix(),
            "taxonomy_file": staged_taxonomy.relative_to(staging).as_posix(),
            "template_lookup": lookup.relative_to(staging).as_posix(),
            "template_count": template_count,
            "unique_template_accession_count": accession_count,
            "taxonomy_source_audit": taxonomy_audit.relative_to(staging).as_posix(),
            "taxonomy_source_counts": taxonomy_counts,
            "components": component_metadata,
            "selection_policy": (
                "SpeciesFinder nominates a bounded candidate set; exact KMA templates "
                "are extracted and must pass BactiPipe ANI/aligned-fraction gates."
            ),
        }
        (staging / "database_identity.json").write_text(
            json.dumps(identity, indent=2, sort_keys=True) + "\n"
        )
        (staging / "BUILD_SUCCESS").write_text(
            f"{built_at}\t{version}\t{status}\n"
        )
        os.replace(staging, final)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return final


def validate_speciesfinder_database(
    path: str | Path, *, allow_candidate: bool = False
) -> SpeciesFinderDatabase:
    root = Path(path).resolve()
    identity_path = root / "database_identity.json"
    success = root / "BUILD_SUCCESS"
    if not identity_path.is_file() or not success.is_file():
        raise ValueError(f"SpeciesFinder database is incomplete: {root}")
    try:
        identity = json.loads(identity_path.read_text())
    except (OSError, ValueError, TypeError) as exc:
        raise ValueError(f"SpeciesFinder database identity is invalid: {identity_path}") from exc
    if identity.get("database_kind") != DATABASE_KIND:
        raise ValueError(f"Not a SpeciesFinder KMA database: {root}")
    if not allow_candidate and not (
        identity.get("status") == "production" and identity.get("active") is True
    ):
        raise ValueError(
            f"Database {identity.get('version', 'unknown')} is not an active production database"
        )
    prefix = root / str(identity.get("database_prefix") or "")
    taxonomy = root / str(identity.get("taxonomy_file") or "")
    lookup = root / str(identity.get("template_lookup") or "")
    _validate_source(prefix, taxonomy)
    if not lookup.is_file() or lookup.stat().st_size == 0:
        raise ValueError(f"SpeciesFinder template lookup is missing: {lookup}")
    components = identity.get("components")
    if not isinstance(components, dict):
        raise ValueError("SpeciesFinder database identity has no component inventory")
    for relative, metadata in components.items():
        component = root / str(relative)
        if not component.is_file() or not isinstance(metadata, dict):
            raise ValueError(f"SpeciesFinder database component is missing: {component}")
        if component.stat().st_size != int(metadata.get("size_bytes", -1)):
            raise ValueError(f"SpeciesFinder database component size changed: {component}")
    return SpeciesFinderDatabase(root, prefix, taxonomy, lookup, identity)


def extract_candidate_references(
    *,
    database: SpeciesFinderDatabase,
    hits: Sequence[SpeciesFinderHit],
    output_directory: str | Path,
    maximum_species: int = 20,
    references_per_species: int = 3,
    kma_executable: str = "kma",
) -> tuple[ExtractedReference, ...]:
    """Extract only nominated KMA templates and publish a local ANI manifest."""

    selected_hits = select_candidate_hits(
        hits,
        maximum_species=maximum_species,
        references_per_species=references_per_species,
    )
    output = Path(output_directory).resolve()
    references = output / "references"
    references.mkdir(parents=True, exist_ok=True)
    if not selected_hits:
        _write_reference_manifest(output / "reference_manifest.tsv", ())
        return ()

    resolved = _resolve_template_ids(database.template_lookup, selected_hits)
    ids = ",".join(str(item[1]) for item in resolved)
    combined = output / "extracted_templates.fna"
    temporary = combined.with_suffix(combined.suffix + ".tmp")
    with temporary.open("w") as stdout:
        completed = subprocess.run(
            [
                kma_executable,
                "seq2fasta",
                "-t_db",
                str(database.database_prefix),
                "-seqs",
                ids,
            ],
            stdout=stdout,
            stderr=subprocess.PIPE,
            text=True,
        )
    (output / "kma_seq2fasta.stderr.log").write_text(completed.stderr or "")
    if completed.returncode != 0:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(
            f"KMA reference extraction failed with exit {completed.returncode}; see "
            f"{output / 'kma_seq2fasta.stderr.log'}"
        )
    os.replace(temporary, combined)
    fasta_records = _read_fasta_records(combined)
    ordered = sorted(resolved, key=lambda item: item[1])
    if len(fasta_records) != len(ordered):
        raise RuntimeError(
            "KMA reference extraction returned a different number of templates "
            f"({len(fasta_records)}) than requested ({len(ordered)})"
        )

    reporting_groups = ReportingGroupCatalog.from_package()
    extracted: list[ExtractedReference] = []
    for (hit, template_id, _header), (fasta_header, sequence) in zip(
        ordered, fasta_records
    ):
        observed_accession = fasta_header.split(" ", 1)[0]
        if observed_accession != hit.template_accession:
            raise RuntimeError(
                "KMA returned an unexpected template: requested "
                f"{hit.template_accession}, received {observed_accession}"
            )
        unique_accession = f"{hit.template_accession}__kma{template_id}"
        destination = references / f"{unique_accession}.fna"
        destination.write_text(f">{fasta_header}\n{sequence}\n")
        extracted.append(
            ExtractedReference(
                accession=unique_accession,
                species_taxid=hit.species_taxid,
                species_name=hit.species_name,
                fasta_path=destination,
                source_template_accession=hit.template_accession,
                template_id=template_id,
                reporting_group=reporting_groups.display_name_for_taxid(
                    hit.species_taxid
                ),
            )
        )
    _write_reference_manifest(output / "reference_manifest.tsv", extracted)
    return tuple(extracted)


def audit_speciesfinder_taxonomy(
    taxonomy_file: str | Path, output_path: str | Path
) -> dict[str, int]:
    """Audit every official template's callable taxonomy/source disposition."""

    source = Path(taxonomy_file)
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    counts: dict[str, int] = {}
    with source.open(newline="", errors="replace") as input_handle, temporary.open(
        "w", newline=""
    ) as output_handle:
        reader = csv.reader(input_handle, delimiter="\t")
        fieldnames = (
            "template_accession",
            "assembly_accession",
            "organism_taxid",
            "species_taxid",
            "species_name",
            "disposition",
            "description",
        )
        writer = csv.DictWriter(output_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for line_number, row in enumerate(reader, start=1):
            if len(row) != 6:
                raise ValueError(
                    f"SpeciesFinder taxonomy line {line_number} has {len(row)} "
                    "columns; expected 6"
                )
            header, assembly, organism_taxid, _taxonomy, species_taxid, species = (
                value.strip() for value in row
            )
            accession, _, description = header.partition(" ")
            if is_excluded_reference_source(description):
                disposition = "excluded_non_isolate_source"
            elif not species or species.casefold() in {"na", "n/a", "unknown"}:
                disposition = "excluded_unresolved_species_name"
            elif not species_taxid.isdigit():
                disposition = "excluded_unresolved_species_taxid"
            elif not is_reportable_species_name(species):
                disposition = "excluded_unresolved_species_label"
            else:
                disposition = "eligible_nomination_template"
            counts[disposition] = counts.get(disposition, 0) + 1
            writer.writerow(
                {
                    "template_accession": accession,
                    "assembly_accession": assembly,
                    "organism_taxid": organism_taxid,
                    "species_taxid": species_taxid,
                    "species_name": species,
                    "disposition": disposition,
                    "description": description,
                }
            )
    os.replace(temporary, destination)
    return counts


def _resolve_template_ids(
    lookup: Path, hits: Sequence[SpeciesFinderHit]
) -> list[tuple[SpeciesFinderHit, int, str]]:
    resolved: list[tuple[SpeciesFinderHit, int, str]] = []
    connection = sqlite3.connect(lookup)
    try:
        for hit in hits:
            rows = connection.execute(
                "SELECT template_id, header FROM templates WHERE accession = ? "
                "ORDER BY template_id",
                (hit.template_accession,),
            ).fetchall()
            if not rows:
                raise ValueError(
                    "SpeciesFinder nominated an accession absent from its KMA name "
                    f"index: {hit.template_accession}"
                )
            template_id, header = rows[0]
            resolved.append((hit, int(template_id), str(header)))
    finally:
        connection.close()
    return resolved


def _write_reference_manifest(
    path: Path, records: Sequence[ExtractedReference]
) -> None:
    canonical_seen: set[str] = set()
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t")
        writer.writeheader()
        for record in records:
            canonical = record.species_taxid not in canonical_seen
            canonical_seen.add(record.species_taxid)
            writer.writerow(
                {
                    "accession": record.accession,
                    "species_taxid": record.species_taxid,
                    "species_name": record.species_name,
                    "fasta_path": record.fasta_path.relative_to(path.parent).as_posix(),
                    "source": "SpeciesFinder KMA template",
                    "reference_role": "representative",
                    "assembly_level": "SpeciesFinder database template",
                    "sequence_origin": "database_reference",
                    "type_material": "no",
                    "canonical_for_tree": "yes" if canonical else "no",
                    "callable_for_identification": "yes",
                    "include_in_kraken": "no",
                    "reporting_group": record.reporting_group,
                    "ncbi_taxonomy_name": _ncbi_species_parent_name(record),
                    "gtdb_taxonomy": "",
                    "sha256": _sha256(record.fasta_path),
                }
            )


def _build_template_lookup(name_file: Path, destination: Path) -> tuple[int, int]:
    connection = sqlite3.connect(destination)
    try:
        connection.execute(
            "CREATE TABLE templates (template_id INTEGER PRIMARY KEY, accession TEXT NOT NULL, header TEXT NOT NULL)"
        )
        connection.execute("CREATE INDEX templates_accession_idx ON templates(accession)")
        rows = []
        accessions: set[str] = set()
        with name_file.open(errors="replace") as handle:
            for template_id, line in enumerate(handle, start=1):
                header = line.rstrip("\r\n")
                accession = header.split(" ", 1)[0]
                if not accession:
                    raise ValueError(f"Blank KMA template name at line {template_id}")
                rows.append((template_id, accession, header))
                accessions.add(accession)
                if len(rows) >= 10_000:
                    connection.executemany(
                        "INSERT INTO templates VALUES (?, ?, ?)", rows
                    )
                    rows.clear()
        if rows:
            connection.executemany("INSERT INTO templates VALUES (?, ?, ?)", rows)
        connection.commit()
        count = int(connection.execute("SELECT COUNT(*) FROM templates").fetchone()[0])
    finally:
        connection.close()
    if count == 0:
        raise ValueError("SpeciesFinder KMA name file contains no templates")
    return count, len(accessions)


def _ncbi_species_parent_name(record: ExtractedReference) -> str:
    if record.species_taxid in {"40041", "148942"}:
        return "Streptococcus equi"
    return record.species_name


def _read_fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    sequence: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(sequence)))
                header = line[1:].strip()
                sequence = []
            else:
                sequence.append(line.strip())
    if header:
        records.append((header, "".join(sequence)))
    if any(not header or not sequence for header, sequence in records):
        raise ValueError(f"KMA produced an invalid FASTA: {path}")
    return records


def _validate_source(prefix: Path, taxonomy: Path) -> None:
    for suffix in REQUIRED_KMA_SUFFIXES:
        component = Path(f"{prefix}{suffix}")
        if not component.is_file() or component.stat().st_size == 0:
            raise ValueError(f"SpeciesFinder database component is missing: {component}")
    if not taxonomy.is_file() or taxonomy.stat().st_size == 0:
        raise ValueError(f"SpeciesFinder taxonomy file is missing: {taxonomy}")


def _copy_or_link(source: Path, destination: Path) -> None:
    try:
        os.link(source, destination)
    except OSError:
        shutil.copy2(source, destination)


def _read_sha256(path: Path) -> str:
    if not path.is_file():
        raise ValueError(f"Source archive checksum file is missing: {path}")
    value = path.read_text().split()[0].lower()
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
        raise ValueError(f"Invalid SHA-256 file: {path}")
    return value


def _speciesfinder_version() -> str:
    try:
        return package_version("SpeciesFinder")
    except PackageNotFoundError:
        return "4.0.4"


def _tool_version(command: Sequence[str]) -> str:
    try:
        completed = subprocess.run(
            list(command), text=True, capture_output=True, timeout=30
        )
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    text = " ".join((completed.stdout or "", completed.stderr or "")).strip()
    return text.splitlines()[0] if text else "unknown"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
