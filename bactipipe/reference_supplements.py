"""Audited recovery of controlled species absent from the GTDB proposal."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

from bactipipe.reference_curation import read_controlled_species
from bactipipe.reference_inventory import AUDIT_COLUMNS


def augment_inventory_with_gap_candidates(
    inventory_tsv: str | Path,
    candidate_assessment_tsv: str | Path,
    controlled_catalog_path: str | Path,
    controlled_gaps_tsv: str | Path,
    output_root: str | Path,
) -> dict[str, object]:
    """Append only quality-passing, exact-taxid candidates for controlled gaps."""

    inventory = _read_tsv(inventory_tsv)
    candidates = _read_tsv(candidate_assessment_tsv)
    gaps = _read_tsv(controlled_gaps_tsv)
    controlled = {record.taxid: record for record in read_controlled_species(controlled_catalog_path)}
    existing_accessions = {row["accession"] for row in inventory}
    existing_taxids = {row["species_taxid"] for row in inventory}
    gap_taxids = {row["species_taxid"] for row in gaps}
    supplements = []
    for candidate in candidates:
        if candidate.get("decision") != "selected_canonical":
            continue
        taxid = candidate.get("species_taxid", "")
        accession = candidate.get("accession", "")
        record = controlled.get(taxid)
        if record is None or taxid not in gap_taxids:
            raise ValueError(
                f"Selected gap candidate {accession} does not match a controlled gap taxid"
            )
        if taxid in existing_taxids:
            raise ValueError(f"Gap candidate taxid {taxid} already exists in the inventory")
        if accession in existing_accessions:
            raise ValueError(f"Gap candidate accession {accession} already exists")
        supplements.append(
            {
                "proposed_accession": accession,
                "accession": accession,
                "accession_resolution_status": "controlled_gap_recovery",
                "species_taxid": taxid,
                "source_species_taxid": taxid,
                "source_taxid": taxid,
                "ncbi_species_name": candidate.get("species_name", ""),
                "gtdb_species_name": "",
                "ncbi_organism_name": candidate.get("submitted_organism", ""),
                "ncbi_genus": _genus(candidate.get("species_name", "")),
                "gtdb_genus": "",
                "source": candidate.get("source", ""),
                "assembly_level": candidate.get("assembly_level", ""),
                "refseq_category": candidate.get("refseq_category", ""),
                "type_material": candidate.get("type_material", "False"),
                "gtdb_representative": "False",
                "sequence_origin": candidate.get("sequence_origin", "isolate"),
                "completeness": candidate.get("completeness", ""),
                "contamination": candidate.get("contamination", ""),
                "contigs": candidate.get("contigs", ""),
                "n50": candidate.get("n50", ""),
                "total_length": candidate.get("total_length", ""),
                "controlled_target": "True",
                "controlled_catalog_taxid": taxid,
                "controlled_catalog_name": record.species_name,
                "controlled_match": "controlled_gap_taxid",
                "reference_role": "canonical",
                "quality_score": candidate.get("quality_score", ""),
                "reason": "Controlled-gap recovery; passed current NCBI quality and taxonomy gates.",
                "current_organism_taxid": taxid,
                "current_organism_name": candidate.get("submitted_organism", ""),
                "current_source": candidate.get("source", ""),
                "current_assembly_level": candidate.get("assembly_level", ""),
                "current_refseq_category": candidate.get("refseq_category", ""),
                "current_type_material": candidate.get("type_material", "False"),
                "current_total_length": candidate.get("total_length", ""),
                "current_checkm_species_taxid": taxid,
                "current_submitted_species": candidate.get("species_name", ""),
                "current_taxonomy_check_status": candidate.get("taxonomy_check_status", ""),
                "metadata_eligible": "True",
                "metadata_exclusion_reason": "",
                "taxonomy_review_required": "False",
                "taxonomy_review_reason": "",
            }
        )
        existing_accessions.add(accession)
        existing_taxids.add(taxid)

    combined = sorted(
        (*inventory, *supplements),
        key=lambda row: (row.get("ncbi_species_name", ""), row.get("reference_role", ""), row["accession"]),
    )
    recovered_taxids = {row["species_taxid"] for row in supplements}
    unresolved = [row for row in gaps if row["species_taxid"] not in recovered_taxids]
    output = Path(output_root)
    output.mkdir(parents=True, exist_ok=True)
    _atomic_tsv(output / "reference_inventory_augmented.tsv", AUDIT_COLUMNS, combined)
    _atomic_tsv(
        output / "accepted_gap_supplements.tsv", AUDIT_COLUMNS, supplements
    )
    _atomic_tsv(
        output / "unresolved_controlled_species.tsv",
        ["stage", "species_taxid", "species_name", "reason"],
        unresolved,
    )
    summary: dict[str, object] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "base_reference_count": len(inventory),
        "supplement_reference_count": len(supplements),
        "final_reference_count": len(combined),
        "final_species_count": len({row["species_taxid"] for row in combined}),
        "controlled_gap_count_before_recovery": len(gaps),
        "controlled_gap_count_after_recovery": len(unresolved),
        "recovered_controlled_taxids": sorted(recovered_taxids, key=int),
        "input_sha256": {
            "inventory": _sha256(Path(inventory_tsv)),
            "candidate_assessment": _sha256(Path(candidate_assessment_tsv)),
            "controlled_catalog": _sha256(Path(controlled_catalog_path)),
            "controlled_gaps": _sha256(Path(controlled_gaps_tsv)),
        },
    }
    _atomic_text(
        output / "augmentation_summary.json",
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
    )
    _atomic_text(
        output / "AUGMENTATION_SUCCESS",
        f"{summary['created_at']}\taugmented_inventory_ready\n",
    )
    return summary


def _genus(name: str) -> str:
    words = name.replace("[", "").replace("]", "").split()
    return words[1] if words and words[0].casefold() == "candidatus" else words[0] if words else ""


def _read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _atomic_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
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
