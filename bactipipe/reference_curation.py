"""Auditable first-pass selection from NCBI Datasets assembly reports."""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import re
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


NON_ISOLATE_PATTERN = re.compile(
    r"(?:^|[^a-z])(mag|metagenome[- ]assembled|metagenomic|sag|single[- ]cell)(?:[^a-z]|$)",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class CurationPolicy:
    minimum_completeness: float = 90.0
    maximum_contamination: float = 5.0
    maximum_contigs: int = 500
    minimum_n50: int = 20_000
    minimum_length: int = 100_000

    def __post_init__(self) -> None:
        if not 0 <= self.minimum_completeness <= 100:
            raise ValueError("minimum_completeness must be between 0 and 100")
        if not 0 <= self.maximum_contamination <= 100:
            raise ValueError("maximum_contamination must be between 0 and 100")
        for name in ("maximum_contigs", "minimum_n50", "minimum_length"):
            if getattr(self, name) < 1:
                raise ValueError(f"{name} must be at least 1")


@dataclass(frozen=True)
class AssemblyCandidate:
    accession: str
    species_taxid: str
    species_name: str
    submitted_organism: str
    source: str
    assembly_level: str
    refseq_category: str
    type_material: bool
    sequence_origin: str
    completeness: float | None
    contamination: float | None
    contigs: int | None
    n50: int | None
    total_length: int | None
    taxonomy_check_status: str
    quality_score: float
    decision: str
    reason: str


@dataclass(frozen=True)
class ControlledSpecies:
    taxid: str
    expected_genome_size: int | None
    legacy_sequence_accession: str
    species_name: str

    @property
    def genus(self) -> str:
        return _name_genus(self.species_name)


@dataclass(frozen=True)
class InventoryReference:
    accession: str
    species_taxid: str
    source_species_taxid: str
    source_taxid: str
    ncbi_species_name: str
    gtdb_species_name: str
    ncbi_organism_name: str
    ncbi_genus: str
    gtdb_genus: str
    source: str
    assembly_level: str
    refseq_category: str
    type_material: bool
    gtdb_representative: bool
    sequence_origin: str
    completeness: float | None
    contamination: float | None
    contigs: int | None
    n50: int | None
    total_length: int | None
    controlled_target: bool
    controlled_catalog_taxid: str
    controlled_catalog_name: str
    controlled_match: str
    reference_role: str
    quality_score: float
    reason: str


INVENTORY_COLUMNS = list(InventoryReference.__dataclass_fields__)
MISSING_CONTROLLED_COLUMNS = [
    "taxid",
    "expected_genome_size",
    "legacy_sequence_accession",
    "species_name",
    "metadata_rows",
    "exclusion_reasons",
]
CONTROLLED_ALIAS_COLUMNS = [
    "current_taxid",
    "catalog_taxids",
    "catalog_names",
    "relationship",
]


def read_controlled_species(path: str | Path) -> tuple[ControlledSpecies, ...]:
    """Read the existing controlled bacterial catalog without changing it."""

    records: list[ControlledSpecies] = []
    seen_taxids: set[str] = set()
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise ValueError(
                    f"Controlled catalog line {line_number} must contain four tab-separated fields"
                )
            taxid, size_text, accession, species_name = (item.strip() for item in fields)
            if not taxid.isdigit() or not species_name:
                raise ValueError(f"Invalid controlled catalog line {line_number}")
            if taxid in seen_taxids:
                raise ValueError(f"Duplicate controlled taxid {taxid} on line {line_number}")
            seen_taxids.add(taxid)
            records.append(
                ControlledSpecies(
                    taxid=taxid,
                    expected_genome_size=_integer(size_text),
                    legacy_sequence_accession=accession,
                    species_name=" ".join(species_name.split()),
                )
            )
    if not records:
        raise ValueError("Controlled species catalog is empty")
    return tuple(records)


def plan_gtdb_reference_inventory(
    metadata_path: str | Path,
    controlled_catalog_path: str | Path,
    *,
    policy: CurationPolicy = CurationPolicy(),
    maximum_diversity_references: int = 3,
    ncbi_taxonomy_dir: str | Path | None = None,
) -> tuple[tuple[InventoryReference, ...], dict[str, object]]:
    """Select a bounded veterinary-plus-near-neighbor reference inventory.

    All assemblies may compete to represent a controlled species. For the
    broader same-genus backbone, only GTDB species representatives may enter.
    This avoids turning a common species with hundreds of thousands of NCBI
    assemblies into an unbounded download request.
    """

    if maximum_diversity_references < 0:
        raise ValueError("maximum_diversity_references cannot be negative")
    controlled = read_controlled_species(controlled_catalog_path)
    merged_taxids = read_merged_taxids(ncbi_taxonomy_dir) if ncbi_taxonomy_dir else {}
    controlled_by_current_taxid: dict[str, list[ControlledSpecies]] = defaultdict(list)
    for record in controlled:
        controlled_by_current_taxid[resolve_merged_taxid(record.taxid, merged_taxids)].append(
            record
        )
    controlled_genera = {record.genus for record in controlled if record.genus}
    candidates: dict[str, list[InventoryReference]] = defaultdict(list)
    counts: Counter[str] = Counter()
    controlled_row_counts: Counter[str] = Counter()
    controlled_exclusions: dict[str, Counter[str]] = defaultdict(Counter)

    with _open_text(metadata_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "accession",
            "gtdb_representative",
            "gtdb_taxonomy",
            "ncbi_taxid",
            "ncbi_species_taxid",
            "ncbi_taxonomy",
            "ncbi_organism_name",
        }
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(
                "GTDB metadata is missing required columns: " + ", ".join(sorted(missing))
            )
        for row in reader:
            counts["metadata_rows"] += 1
            source_species_taxid = str(row.get("ncbi_species_taxid") or "").strip()
            source_taxid = str(row.get("ncbi_taxid") or "").strip()
            current_species_taxid = resolve_merged_taxid(
                source_species_taxid, merged_taxids
            )
            current_source_taxid = resolve_merged_taxid(source_taxid, merged_taxids)
            ncbi_species = _taxonomy_rank(row.get("ncbi_taxonomy"), "s")
            gtdb_species = _taxonomy_rank(row.get("gtdb_taxonomy"), "s")
            ncbi_genus = _taxonomy_rank(row.get("ncbi_taxonomy"), "g")
            gtdb_genus = _taxonomy_rank(row.get("gtdb_taxonomy"), "g")
            normalized_ncbi_species = _normalize_species_name(ncbi_species)
            controlled_match = ""
            specific_matches = controlled_by_current_taxid.get(current_source_taxid, [])
            species_matches = controlled_by_current_taxid.get(current_species_taxid, [])
            controlled_matches = specific_matches or species_matches
            controlled_catalog_taxid = ";".join(
                record.taxid for record in controlled_matches
            )
            controlled_catalog_name = ";".join(
                record.species_name for record in controlled_matches
            )
            if controlled_matches:
                match_field = "ncbi_taxid" if specific_matches else "ncbi_species_taxid"
                merged = any(
                    resolve_merged_taxid(record.taxid, merged_taxids) != record.taxid
                    for record in controlled_matches
                )
                controlled_match = (
                    f"merged_{match_field}" if merged else match_field
                )
            controlled_target = bool(controlled_match)
            for controlled_record in controlled_matches:
                controlled_row_counts[controlled_record.taxid] += 1
            is_representative = _truthy(row.get("gtdb_representative"))
            is_near_neighbor = is_representative and bool(
                {ncbi_genus.casefold(), gtdb_genus.casefold()}.intersection(
                    controlled_genera
                )
            )
            if not controlled_target and not is_near_neighbor:
                continue
            counts["relevant_rows"] += 1
            if not ncbi_species:
                if controlled_matches:
                    ncbi_species = controlled_matches[0].species_name
                    normalized_ncbi_species = _normalize_species_name(ncbi_species)
                else:
                    counts["excluded:unresolved_ncbi_species_name"] += 1
                    continue
            exclusion = _gtdb_exclusion_reason(row, policy=policy)
            if exclusion:
                counts[f"excluded:{exclusion}"] += 1
                for controlled_record in controlled_matches:
                    controlled_exclusions[controlled_record.taxid][exclusion] += 1
                continue
            accession = _normalize_assembly_accession(row.get("accession"))
            if not accession:
                counts["excluded:invalid_assembly_accession"] += 1
                for controlled_record in controlled_matches:
                    controlled_exclusions[controlled_record.taxid][
                        "invalid_assembly_accession"
                    ] += 1
                continue
            specific_controlled_taxid = ""
            if specific_matches and current_source_taxid != current_species_taxid:
                specific_controlled_taxid = current_source_taxid
            reference_taxid = specific_controlled_taxid or current_species_taxid
            species_key = reference_taxid if reference_taxid.isdigit() else normalized_ncbi_species
            if not species_key:
                counts["excluded:unresolved_ncbi_species"] += 1
                continue
            completeness = _preferred_number(row, "checkm2_completeness", "checkm_completeness")
            contamination = _preferred_number(row, "checkm2_contamination", "checkm_contamination")
            type_material = _has_type_material(row.get("ncbi_type_material_designation"))
            reference = InventoryReference(
                accession=accession,
                species_taxid=reference_taxid,
                source_species_taxid=source_species_taxid,
                source_taxid=source_taxid,
                ncbi_species_name=ncbi_species,
                gtdb_species_name=gtdb_species,
                ncbi_organism_name=" ".join(str(row.get("ncbi_organism_name") or "").split()),
                ncbi_genus=ncbi_genus,
                gtdb_genus=gtdb_genus,
                source="RefSeq" if accession.startswith("GCF_") else "GenBank",
                assembly_level=str(row.get("ncbi_assembly_level") or "").strip(),
                refseq_category=str(row.get("ncbi_refseq_category") or "").strip(),
                type_material=type_material,
                gtdb_representative=is_representative,
                sequence_origin="isolate",
                completeness=completeness,
                contamination=contamination,
                contigs=_preferred_integer(row, "ncbi_contig_count", "contig_count"),
                n50=_preferred_integer(row, "ncbi_contig_n50", "n50_contigs"),
                total_length=_preferred_integer(row, "ncbi_total_length", "genome_size"),
                controlled_target=controlled_target,
                controlled_catalog_taxid=controlled_catalog_taxid,
                controlled_catalog_name=controlled_catalog_name,
                controlled_match=controlled_match,
                reference_role="candidate",
                quality_score=round(_gtdb_quality_score(row), 6),
                reason=(
                    "Controlled veterinary species candidate."
                    if controlled_target
                    else "GTDB representative from a controlled veterinary genus."
                ),
            )
            group = candidates[species_key]
            group.append(reference)
            # Common controlled species may contribute enormous numbers of
            # assemblies. Retain only a bounded high-quality pool while always
            # preserving GTDB representatives as diversity candidates.
            if len(group) > 24:
                representatives = [item for item in group if item.gtdb_representative]
                others = sorted(
                    (item for item in group if not item.gtdb_representative),
                    key=_inventory_rank,
                    reverse=True,
                )[:8]
                candidates[species_key] = representatives + others

    selected: list[InventoryReference] = []
    for group in candidates.values():
        ordered = sorted(group, key=_inventory_rank, reverse=True)
        canonical = ordered[0]
        selected.append(
            InventoryReference(
                **{
                    **asdict(canonical),
                    "reference_role": "canonical",
                    "reason": canonical.reason + " Highest-ranked eligible assembly.",
                }
            )
        )
        diversity = [
            item
            for item in ordered[1:]
            if item.gtdb_representative and item.accession != canonical.accession
        ][:maximum_diversity_references]
        for item in diversity:
            selected.append(
                InventoryReference(
                    **{
                        **asdict(item),
                        "reference_role": "diversity",
                        "reason": item.reason + " Retained GTDB-cluster diversity representative.",
                    }
                )
            )

    unique = {item.accession: item for item in selected}
    final = tuple(
        sorted(
            unique.values(),
            key=lambda item: (item.ncbi_species_name, item.reference_role, item.accession),
        )
    )
    counts["selected_references"] = len(final)
    counts["selected_species"] = len({item.species_taxid or item.ncbi_species_name for item in final})
    counts["controlled_catalog_species"] = len(controlled)
    represented_controlled_taxids = {
        taxid
        for item in final
        for taxid in item.controlled_catalog_taxid.split(";")
        if taxid
    }
    missing_controlled = []
    for record in controlled:
        if record.taxid in represented_controlled_taxids:
            continue
        exclusions = controlled_exclusions.get(record.taxid, Counter())
        missing_controlled.append(
            {
                **asdict(record),
                "metadata_rows": controlled_row_counts[record.taxid],
                "exclusion_reasons": ";".join(
                    f"{reason}={count}" for reason, count in sorted(exclusions.items())
                )
                or "no_matching_gtdb_metadata",
            }
        )
    counts["controlled_species_represented"] = len(represented_controlled_taxids)
    counts["controlled_species_missing"] = len(missing_controlled)
    controlled_aliases = []
    for current_taxid, records in sorted(
        controlled_by_current_taxid.items(), key=lambda item: int(item[0])
    ):
        merged_records = [
            record
            for record in records
            if resolve_merged_taxid(record.taxid, merged_taxids) != record.taxid
        ]
        if len(records) < 2 and not merged_records:
            continue
        controlled_aliases.append(
            {
                "current_taxid": current_taxid,
                "catalog_taxids": ";".join(record.taxid for record in records),
                "catalog_names": ";".join(record.species_name for record in records),
                "relationship": (
                    "multiple_catalog_entries_resolve_to_one_current_taxid"
                    if len(records) > 1
                    else "obsolete_catalog_taxid_resolved"
                ),
            }
        )
    counts["controlled_taxid_alias_group_count"] = len(controlled_aliases)
    summary: dict[str, object] = {
        "schema_version": 1,
        "selection_scope": (
            "All eligible assemblies for controlled species plus eligible GTDB "
            "representatives from the same NCBI or GTDB genera"
        ),
        "maximum_diversity_references_per_species": maximum_diversity_references,
        "ncbi_merged_taxid_count": len(merged_taxids),
        "counts": dict(sorted(counts.items())),
        "missing_controlled_species": missing_controlled,
        "controlled_taxid_aliases": controlled_aliases,
    }
    return final, summary


def write_inventory_plan(
    references: Iterable[InventoryReference],
    *,
    output_tsv: str | Path,
    accessions_path: str | Path,
    summary: dict[str, object],
    summary_path: str | Path,
    missing_controlled_path: str | Path | None = None,
    controlled_aliases_path: str | Path | None = None,
) -> None:
    references = tuple(references)
    _atomic_tsv(output_tsv, INVENTORY_COLUMNS, (asdict(item) for item in references))
    _atomic_text(
        accessions_path,
        "".join(f"{item.accession}\n" for item in references),
    )
    _atomic_text(summary_path, json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if missing_controlled_path is not None:
        missing = summary.get("missing_controlled_species") or []
        if not isinstance(missing, list):
            raise ValueError("Inventory summary has invalid missing_controlled_species")
        _atomic_tsv(
            missing_controlled_path,
            MISSING_CONTROLLED_COLUMNS,
            missing,
        )
    if controlled_aliases_path is not None:
        aliases = summary.get("controlled_taxid_aliases") or []
        if not isinstance(aliases, list):
            raise ValueError("Inventory summary has invalid controlled_taxid_aliases")
        _atomic_tsv(controlled_aliases_path, CONTROLLED_ALIAS_COLUMNS, aliases)


def plan_references(
    reports: Iterable[dict[str, object]],
    *,
    policy: CurationPolicy = CurationPolicy(),
) -> tuple[AssemblyCandidate, ...]:
    assessed = [_assess(record, policy=policy) for record in reports]
    eligible: dict[str, list[AssemblyCandidate]] = {}
    for candidate in assessed:
        if candidate.decision == "eligible":
            eligible.setdefault(candidate.species_taxid, []).append(candidate)
    selected: dict[str, str] = {}
    for taxid, candidates in eligible.items():
        winner = max(candidates, key=lambda item: (item.quality_score, item.accession))
        selected[taxid] = winner.accession
    final = []
    for candidate in assessed:
        if selected.get(candidate.species_taxid) == candidate.accession:
            final.append(
                AssemblyCandidate(
                    **{
                        **asdict(candidate),
                        "decision": "selected_canonical",
                        "reason": "Highest-ranked eligible assembly for this species.",
                    }
                )
            )
        elif candidate.decision == "eligible":
            final.append(
                AssemblyCandidate(
                    **{
                        **asdict(candidate),
                        "decision": "eligible_not_selected",
                        "reason": "A higher-ranked eligible assembly represents this species.",
                    }
                )
            )
        else:
            final.append(candidate)
    return tuple(sorted(final, key=lambda item: (item.species_name, item.accession)))


def read_ncbi_jsonl(path: str | Path) -> Iterable[dict[str, object]]:
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            value = json.loads(line)
            if not isinstance(value, dict):
                raise ValueError(f"NCBI report line {line_number} is not a JSON object")
            yield value


def read_merged_taxids(ncbi_taxonomy_dir: str | Path) -> dict[str, str]:
    """Read NCBI's authoritative obsolete-to-current taxid mapping."""

    source = Path(ncbi_taxonomy_dir) / "merged.dmp"
    if not source.is_file():
        raise ValueError(f"NCBI taxonomy merged.dmp was not found: {source}")
    merged: dict[str, str] = {}
    with source.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = [field.strip() for field in line.split("|")]
            if len(fields) < 2 or not fields[0].isdigit() or not fields[1].isdigit():
                raise ValueError(f"Invalid merged.dmp record on line {line_number}")
            previous, current = fields[:2]
            if previous in merged:
                raise ValueError(f"Duplicate obsolete taxid {previous} in merged.dmp")
            merged[previous] = current
    return merged


def resolve_merged_taxid(taxid: object, merged: dict[str, str]) -> str:
    current = str(taxid or "").strip()
    visited: set[str] = set()
    while current in merged:
        if current in visited:
            raise ValueError(f"Cycle detected in NCBI merged taxids at {current}")
        visited.add(current)
        current = merged[current]
    return current


def write_curation_plan(candidates: Iterable[AssemblyCandidate], path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    rows = [asdict(candidate) for candidate in candidates]
    fields = list(AssemblyCandidate.__dataclass_fields__)
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, destination)


def _assess(record: dict[str, object], *, policy: CurationPolicy) -> AssemblyCandidate:
    assembly_info = _dict(record.get("assembly_info"))
    stats = _dict(record.get("assembly_stats"))
    checkm = _dict(record.get("checkm_info"))
    ani = _dict(record.get("average_nucleotide_identity"))
    organism = _dict(record.get("organism"))
    accession = str(record.get("accession") or "").strip()
    species_taxid = str(checkm.get("checkm_species_tax_id") or "").strip()
    species_name = " ".join(str(ani.get("submitted_species") or "").split())
    submitted = " ".join(
        str(ani.get("submitted_organism") or organism.get("organism_name") or "").split()
    )
    source = str(record.get("source_database") or "").replace("SOURCE_DATABASE_", "")
    assembly_level = str(assembly_info.get("assembly_level") or "")
    refseq_category = str(assembly_info.get("refseq_category") or "")
    type_material = bool(record.get("type_material"))
    completeness = _number(checkm.get("completeness"))
    contamination = _number(checkm.get("contamination"))
    contigs = _integer(stats.get("number_of_contigs"))
    n50 = _integer(stats.get("contig_n50"))
    total_length = _integer(stats.get("total_sequence_length"))
    taxonomy_status = str(ani.get("taxonomy_check_status") or "")
    origin = "non_isolate" if _looks_non_isolate(record) else "isolate"
    failures = []
    if not accession or not species_taxid.isdigit() or not species_name:
        failures.append("unresolved_species_metadata")
    if assembly_info.get("assembly_status") not in {None, "", "current"}:
        failures.append("assembly_not_current")
    if origin != "isolate":
        failures.append("non_isolate_source")
    if taxonomy_status and taxonomy_status != "OK":
        failures.append("ncbi_taxonomy_check_failed")
    if completeness is None or completeness < policy.minimum_completeness:
        failures.append("low_or_missing_completeness")
    if contamination is None or contamination > policy.maximum_contamination:
        failures.append("high_or_missing_contamination")
    if contigs is None or contigs > policy.maximum_contigs:
        failures.append("excessive_or_missing_contig_count")
    if n50 is None or n50 < policy.minimum_n50:
        failures.append("low_or_missing_n50")
    if total_length is None or total_length < policy.minimum_length:
        failures.append("short_or_missing_assembly_length")
    quality_score = _quality_score(
        source=source,
        refseq_category=refseq_category,
        assembly_level=assembly_level,
        type_material=type_material,
        completeness=completeness,
        contamination=contamination,
        contigs=contigs,
        n50=n50,
    )
    return AssemblyCandidate(
        accession=accession,
        species_taxid=species_taxid,
        species_name=species_name,
        submitted_organism=submitted,
        source=source,
        assembly_level=assembly_level,
        refseq_category=refseq_category,
        type_material=type_material,
        sequence_origin=origin,
        completeness=completeness,
        contamination=contamination,
        contigs=contigs,
        n50=n50,
        total_length=total_length,
        taxonomy_check_status=taxonomy_status,
        quality_score=round(quality_score, 6),
        decision="excluded" if failures else "eligible",
        reason=";".join(failures) if failures else "Passed automated quality gates.",
    )


def _quality_score(
    *,
    source: str,
    refseq_category: str,
    assembly_level: str,
    type_material: bool,
    completeness: float | None,
    contamination: float | None,
    contigs: int | None,
    n50: int | None,
) -> float:
    score = 10_000 if type_material else 0
    if source == "REFSEQ":
        score += 1_000
    category = refseq_category.casefold()
    if category == "reference genome":
        score += 5_000
    elif category == "representative genome":
        score += 2_500
    score += {
        "complete genome": 500,
        "chromosome": 400,
        "scaffold": 200,
        "contig": 100,
    }.get(assembly_level.casefold(), 0)
    score += completeness or 0
    score -= 5 * (contamination or 0)
    score -= (contigs or 0) / 100
    if n50 and n50 > 0:
        score += math.log10(n50)
    return score


def _looks_non_isolate(record: dict[str, object]) -> bool:
    text = json.dumps(record, sort_keys=True)
    return bool(NON_ISOLATE_PATTERN.search(text))


def _dict(value: object) -> dict[str, object]:
    return value if isinstance(value, dict) else {}


def _number(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _integer(value: object) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _open_text(path: str | Path):
    source = Path(path)
    if source.suffix == ".gz":
        return gzip.open(source, "rt", newline="")
    return source.open(newline="")


def _taxonomy_rank(value: object, prefix: str) -> str:
    for part in str(value or "").split(";"):
        part = part.strip()
        if part.startswith(f"{prefix}__"):
            return " ".join(part[3:].split())
    return ""


def _normalize_species_name(value: object) -> str:
    text = " ".join(str(value or "").replace("[", "").replace("]", "").split())
    return text.casefold()


def _name_genus(value: object) -> str:
    words = str(value or "").replace("[", "").replace("]", "").split()
    if not words:
        return ""
    if words[0].casefold() == "candidatus" and len(words) > 1:
        return words[1].casefold()
    return words[0].casefold()


def _truthy(value: object) -> bool:
    return str(value or "").strip().casefold() in {"1", "true", "t", "yes", "y"}


def _normalize_assembly_accession(value: object) -> str:
    accession = str(value or "").strip()
    if accession.startswith(("RS_", "GB_")):
        accession = accession[3:]
    return accession if re.fullmatch(r"GC[AF]_\d+\.\d+", accession) else ""


def _preferred_number(row: dict[str, str], primary: str, fallback: str) -> float | None:
    return _number(row.get(primary)) if _number(row.get(primary)) is not None else _number(row.get(fallback))


def _preferred_integer(row: dict[str, str], primary: str, fallback: str) -> int | None:
    preferred = _integer(row.get(primary))
    return preferred if preferred is not None else _integer(row.get(fallback))


def _has_type_material(value: object) -> bool:
    normalized = str(value or "").strip().casefold()
    return bool(normalized and normalized not in {"none", "na", "n/a", "not available"})


def _gtdb_exclusion_reason(row: dict[str, str], *, policy: CurationPolicy) -> str:
    category = str(row.get("ncbi_genome_category") or "").strip().casefold()
    if category not in {"", "none", "na", "n/a"}:
        return "non_isolate_source"
    representation = str(row.get("ncbi_genome_representation") or "").strip().casefold()
    if representation and representation not in {"full", "complete"}:
        return "partial_genome_representation"
    completeness = _preferred_number(row, "checkm2_completeness", "checkm_completeness")
    contamination = _preferred_number(row, "checkm2_contamination", "checkm_contamination")
    contigs = _preferred_integer(row, "ncbi_contig_count", "contig_count")
    n50 = _preferred_integer(row, "ncbi_contig_n50", "n50_contigs")
    total_length = _preferred_integer(row, "ncbi_total_length", "genome_size")
    if completeness is None or completeness < policy.minimum_completeness:
        return "low_or_missing_completeness"
    if contamination is None or contamination > policy.maximum_contamination:
        return "high_or_missing_contamination"
    if contigs is None or contigs > policy.maximum_contigs:
        return "excessive_or_missing_contig_count"
    if n50 is None or n50 < policy.minimum_n50:
        return "low_or_missing_n50"
    if total_length is None or total_length < policy.minimum_length:
        return "short_or_missing_assembly_length"
    return ""


def _gtdb_quality_score(row: dict[str, str]) -> float:
    accession = _normalize_assembly_accession(row.get("accession"))
    completeness = _preferred_number(row, "checkm2_completeness", "checkm_completeness")
    contamination = _preferred_number(row, "checkm2_contamination", "checkm_contamination")
    contigs = _preferred_integer(row, "ncbi_contig_count", "contig_count")
    n50 = _preferred_integer(row, "ncbi_contig_n50", "n50_contigs")
    return _quality_score(
        source="REFSEQ" if accession.startswith("GCF_") else "GENBANK",
        refseq_category=str(row.get("ncbi_refseq_category") or ""),
        assembly_level=str(row.get("ncbi_assembly_level") or ""),
        type_material=_has_type_material(row.get("ncbi_type_material_designation")),
        completeness=completeness,
        contamination=contamination,
        contigs=contigs,
        n50=n50,
    ) + (750 if _truthy(row.get("gtdb_representative")) else 0)


def _inventory_rank(item: InventoryReference) -> tuple[float, int, str]:
    return (item.quality_score, int(item.gtdb_representative), item.accession)


def _atomic_tsv(path: str | Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, destination)


def _atomic_text(path: str | Path, content: str) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.write_text(content)
    os.replace(temporary, destination)
