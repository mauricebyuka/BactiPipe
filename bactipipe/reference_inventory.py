"""Reconcile a GTDB-derived proposal against current NCBI assembly records."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from bactipipe.ncbi_reference_metadata import ASSEMBLY_ACCESSION
from bactipipe.reference_curation import read_merged_taxids, resolve_merged_taxid


AUDIT_COLUMNS = [
    "proposed_accession",
    "accession",
    "accession_resolution_status",
    "species_taxid",
    "source_species_taxid",
    "source_taxid",
    "ncbi_species_name",
    "gtdb_species_name",
    "ncbi_organism_name",
    "ncbi_genus",
    "gtdb_genus",
    "source",
    "assembly_level",
    "refseq_category",
    "type_material",
    "gtdb_representative",
    "sequence_origin",
    "completeness",
    "contamination",
    "contigs",
    "n50",
    "total_length",
    "controlled_target",
    "controlled_catalog_taxid",
    "controlled_catalog_name",
    "controlled_match",
    "reference_role",
    "quality_score",
    "reason",
    "current_organism_taxid",
    "current_organism_name",
    "current_source",
    "current_assembly_level",
    "current_refseq_category",
    "current_type_material",
    "current_total_length",
    "current_checkm_species_taxid",
    "current_submitted_species",
    "current_taxonomy_check_status",
    "metadata_eligible",
    "metadata_exclusion_reason",
    "taxonomy_review_required",
    "taxonomy_review_reason",
]


def reconcile_reference_inventory(
    inventory_tsv: str | Path,
    resolution_tsv: str | Path,
    reports_jsonl: str | Path,
    output_root: str | Path,
    *,
    pre_inventory_missing_tsv: str | Path | None = None,
    ncbi_taxonomy_dir: str | Path | None = None,
) -> dict[str, object]:
    inventory = _read_tsv(inventory_tsv)
    resolution_rows = _read_tsv(resolution_tsv)
    reports = _read_reports(reports_jsonl)
    resolutions = _read_resolutions(resolution_rows)
    merged_taxids = read_merged_taxids(ncbi_taxonomy_dir) if ncbi_taxonomy_dir else {}
    proposed_accessions = {row.get("accession", "") for row in inventory}
    if set(resolutions) != proposed_accessions:
        missing = proposed_accessions.difference(resolutions)
        extra = set(resolutions).difference(proposed_accessions)
        raise ValueError(
            "Resolution/inventory accession mismatch; "
            f"missing={sorted(missing)[:3]}, extra={sorted(extra)[:3]}"
        )

    resolved_by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    missing_rows: list[dict[str, str]] = []
    excluded_metadata_rows: list[dict[str, str]] = []
    returned_seen: set[str] = set()
    review_counts: Counter[str] = Counter()
    for row in inventory:
        proposed = row["accession"]
        resolution = resolutions[proposed]
        returned = resolution["returned_accession"]
        if resolution["status"] == "missing" or not returned:
            missing_rows.append(
                {
                    **row,
                    "proposed_accession": proposed,
                    "missing_reason": "not_returned_by_ncbi_current_query",
                }
            )
            continue
        if returned in returned_seen:
            raise ValueError(f"More than one proposal resolved to {returned}")
        returned_seen.add(returned)
        report = reports.get(returned)
        if report is None:
            raise ValueError(f"Resolved accession lacks an NCBI report: {returned}")
        enriched = _enrich_row(row, resolution, report, merged_taxids=merged_taxids)
        for reason in filter(None, enriched["taxonomy_review_reason"].split(";")):
            review_counts[reason] += 1
        if enriched["metadata_eligible"] != "True":
            excluded_metadata_rows.append(enriched)
            continue
        species_key = row.get("species_taxid") or row.get("ncbi_species_name") or returned
        resolved_by_species[species_key].append(enriched)

    final: list[dict[str, str]] = []
    missing_species: list[dict[str, str]] = []
    all_species = defaultdict(list)
    for row in inventory:
        species_key = row.get("species_taxid") or row.get("ncbi_species_name") or row["accession"]
        all_species[species_key].append(row)
    for species_key, proposed_group in all_species.items():
        resolved_group = resolved_by_species.get(species_key, [])
        if not resolved_group:
            example = proposed_group[0]
            missing_species.append(
                {
                    "species_taxid": example.get("species_taxid", ""),
                    "ncbi_species_name": example.get("ncbi_species_name", ""),
                    "controlled_target": example.get("controlled_target", "False"),
                    "proposed_reference_count": str(len(proposed_group)),
                    "reason": "no_current_metadata_eligible_reference",
                }
            )
            continue
        canonical = next(
            (row for row in resolved_group if row["reference_role"] == "canonical"),
            None,
        )
        if canonical is None:
            canonical = max(
                resolved_group,
                key=lambda row: (float(row.get("quality_score") or 0), row["accession"]),
            )
            canonical["reference_role"] = "canonical"
            canonical["reason"] += " Promoted because the proposed canonical accession was unavailable."
        for row in resolved_group:
            if row is not canonical and row["reference_role"] == "canonical":
                row["reference_role"] = "diversity"
            final.append(row)

    final.sort(key=lambda row: (row["ncbi_species_name"], row["reference_role"], row["accession"]))
    output = Path(output_root)
    output.mkdir(parents=True, exist_ok=True)
    _atomic_tsv(output / "reference_inventory_current.tsv", AUDIT_COLUMNS, final)
    _atomic_text(
        output / "current_assembly_accessions.txt",
        "".join(f"{row['accession']}\n" for row in final),
    )
    missing_columns = [
        "proposed_accession",
        "species_taxid",
        "ncbi_species_name",
        "reference_role",
        "controlled_target",
        "missing_reason",
    ]
    _atomic_tsv(output / "missing_proposed_accessions.tsv", missing_columns, missing_rows)
    _atomic_tsv(
        output / "excluded_current_metadata.tsv",
        AUDIT_COLUMNS,
        excluded_metadata_rows,
    )
    taxonomy_review_rows = [
        row
        for row in (*final, *excluded_metadata_rows)
        if row["taxonomy_review_required"] == "True"
    ]
    taxonomy_review_rows.sort(
        key=lambda row: (row["ncbi_species_name"], row["accession"])
    )
    _atomic_tsv(output / "taxonomy_review.tsv", AUDIT_COLUMNS, taxonomy_review_rows)
    _atomic_tsv(
        output / "missing_species.tsv",
        [
            "species_taxid",
            "ncbi_species_name",
            "controlled_target",
            "proposed_reference_count",
            "reason",
        ],
        missing_species,
    )
    controlled_gaps: list[dict[str, str]] = []
    if pre_inventory_missing_tsv is not None:
        for row in _read_tsv(pre_inventory_missing_tsv):
            controlled_gaps.append(
                {
                    "stage": "gtdb_inventory",
                    "species_taxid": row.get("taxid", ""),
                    "species_name": row.get("species_name", ""),
                    "reason": row.get("exclusion_reasons", ""),
                }
            )
    controlled_gaps.extend(
        {
            "stage": "ncbi_current_resolution",
            "species_taxid": row["species_taxid"],
            "species_name": row["ncbi_species_name"],
            "reason": row["reason"],
        }
        for row in missing_species
        if row["controlled_target"] == "True"
    )
    controlled_gaps.sort(key=lambda row: (row["species_name"], row["species_taxid"]))
    _atomic_tsv(
        output / "controlled_species_gaps.tsv",
        ["stage", "species_taxid", "species_name", "reason"],
        controlled_gaps,
    )
    summary: dict[str, object] = {
        "schema_version": 2,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "proposed_reference_count": len(inventory),
        "current_reference_count": len(final),
        "missing_proposed_accession_count": len(missing_rows),
        "excluded_current_metadata_count": len(excluded_metadata_rows),
        "excluded_current_metadata_reason_counts": dict(
            sorted(
                Counter(
                    row["metadata_exclusion_reason"] for row in excluded_metadata_rows
                ).items()
            )
        ),
        "proposed_species_count": len(all_species),
        "current_species_count": len(resolved_by_species),
        "missing_species_count": len(missing_species),
        "missing_controlled_species_count": sum(
            row["controlled_target"] == "True" for row in missing_species
        ),
        "total_controlled_species_gap_count": len(controlled_gaps),
        "taxonomy_review_reference_count": len(taxonomy_review_rows),
        "taxonomy_review_species_count": len(
            {row["species_taxid"] for row in taxonomy_review_rows}
        ),
        "taxonomy_review_reason_counts": dict(sorted(review_counts.items())),
        "ncbi_merged_taxid_count": len(merged_taxids),
        "input_sha256": {
            "inventory": _sha256(Path(inventory_tsv)),
            "accession_resolution": _sha256(Path(resolution_tsv)),
            "ncbi_reports": _sha256(Path(reports_jsonl)),
            "pre_inventory_missing": (
                _sha256(Path(pre_inventory_missing_tsv))
                if pre_inventory_missing_tsv is not None
                else ""
            ),
            "ncbi_merged_taxids": (
                _sha256(Path(ncbi_taxonomy_dir) / "merged.dmp")
                if ncbi_taxonomy_dir is not None
                else ""
            ),
        },
    }
    _atomic_text(
        output / "reconciliation_summary.json",
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
    )
    _atomic_text(
        output / "RECONCILIATION_SUCCESS",
        f"{summary['created_at']}\tcurrent_inventory_ready_for_review\n",
    )
    return summary


def _read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_reports(path: str | Path) -> dict[str, dict[str, object]]:
    reports: dict[str, dict[str, object]] = {}
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            value = json.loads(line)
            accession = str(value.get("accession") or "") if isinstance(value, dict) else ""
            if not ASSEMBLY_ACCESSION.fullmatch(accession):
                raise ValueError(f"Invalid assembly report on line {line_number}")
            if accession in reports:
                raise ValueError(f"Duplicate assembly report {accession}")
            reports[accession] = value
    return reports


def _read_resolutions(rows: Iterable[dict[str, str]]) -> dict[str, dict[str, str]]:
    resolutions: dict[str, dict[str, str]] = {}
    for row in rows:
        requested = row.get("requested_accession", "")
        if not ASSEMBLY_ACCESSION.fullmatch(requested) or requested in resolutions:
            raise ValueError(f"Invalid or duplicate requested accession: {requested}")
        status = row.get("status", "")
        if status not in {"exact", "current_version_substituted", "missing"}:
            raise ValueError(f"Invalid accession-resolution status: {status}")
        resolutions[requested] = row
    return resolutions


def _enrich_row(
    row: dict[str, str],
    resolution: dict[str, str],
    report: dict[str, object],
    *,
    merged_taxids: dict[str, str],
) -> dict[str, str]:
    assembly_info = _dict(report.get("assembly_info"))
    assembly_stats = _dict(report.get("assembly_stats"))
    organism = _dict(report.get("organism"))
    checkm = _dict(report.get("checkm_info"))
    ani = _dict(report.get("average_nucleotide_identity"))
    proposed = row["accession"]
    returned = resolution["returned_accession"]
    reasons = []
    exclusion_reasons = []
    current_species_taxid = str(checkm.get("checkm_species_tax_id") or "")
    proposed_species_taxid = (
        row.get("source_species_taxid") or row.get("species_taxid", "")
    )
    resolved_current_species_taxid = resolve_merged_taxid(
        current_species_taxid, merged_taxids
    )
    resolved_proposed_species_taxid = resolve_merged_taxid(
        proposed_species_taxid, merged_taxids
    )
    if (
        resolved_current_species_taxid
        and resolved_proposed_species_taxid
        and resolved_current_species_taxid != resolved_proposed_species_taxid
    ):
        reasons.append("current_species_taxid_changed")
    current_organism_taxid = str(organism.get("tax_id") or "")
    proposed_organism_taxid = row.get("source_taxid") or proposed_species_taxid
    if (
        current_organism_taxid
        and proposed_organism_taxid
        and resolve_merged_taxid(current_organism_taxid, merged_taxids)
        != resolve_merged_taxid(proposed_organism_taxid, merged_taxids)
    ):
        reasons.append("current_organism_taxid_changed")
    submitted_species = " ".join(str(ani.get("submitted_species") or "").split())
    species_name_changed = submitted_species and _normalize_name(
        submitted_species
    ) != _normalize_name(row.get("ncbi_species_name", ""))
    if species_name_changed:
        if row.get("controlled_target") != "True" and _is_unnamed_species(
            submitted_species
        ):
            exclusion_reasons.append("unnamed_current_species")
    taxonomy_status = str(ani.get("taxonomy_check_status") or "")
    if taxonomy_status.casefold() == "failed":
        reasons.append("ncbi_taxonomy_check_failed")
    if reasons:
        exclusion_reasons.append("taxonomy_assignment_conflict")
    if assembly_info.get("assembly_status") not in {None, "", "current"}:
        exclusion_reasons.append("assembly_not_current")
    return {
        **row,
        "proposed_accession": proposed,
        "accession": returned,
        "accession_resolution_status": resolution["status"],
        "current_organism_taxid": current_organism_taxid,
        "current_organism_name": " ".join(str(organism.get("organism_name") or "").split()),
        "current_source": str(report.get("source_database") or "").replace(
            "SOURCE_DATABASE_", ""
        ),
        "current_assembly_level": str(assembly_info.get("assembly_level") or ""),
        "current_refseq_category": str(assembly_info.get("refseq_category") or ""),
        "current_type_material": str(bool(report.get("type_material"))),
        "current_total_length": str(assembly_stats.get("total_sequence_length") or ""),
        "current_checkm_species_taxid": current_species_taxid,
        "current_submitted_species": submitted_species,
        "current_taxonomy_check_status": taxonomy_status,
        "metadata_eligible": str(not exclusion_reasons),
        "metadata_exclusion_reason": ";".join(exclusion_reasons),
        "taxonomy_review_required": str(bool(reasons)),
        "taxonomy_review_reason": ";".join(reasons),
    }


def _dict(value: object) -> dict[str, object]:
    return value if isinstance(value, dict) else {}


def _normalize_name(value: object) -> str:
    return " ".join(str(value or "").split()).casefold()


def _is_unnamed_species(value: object) -> bool:
    normalized = f" {_normalize_name(value)} "
    return any(
        marker in normalized
        for marker in (
            " sp. ",
            " bacterium ",
            " uncultured ",
            " endosymbiont ",
            " symbiont ",
        )
    )


def _atomic_tsv(
    path: Path,
    fields: list[str],
    rows: Iterable[dict[str, str]],
) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def _atomic_text(path: Path, content: str) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(content)
    os.replace(temporary, path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
