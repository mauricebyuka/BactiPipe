from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from bactipipe.reference_inventory import reconcile_reference_inventory


class ReferenceInventoryTests(unittest.TestCase):
    def test_reconciliation_promotes_fallback_and_audits_taxonomy_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inventory = root / "inventory.tsv"
            inventory_rows = [
                inventory_row("GCF_0001.1", "1", "Species alpha", "canonical", "True", 100),
                inventory_row("GCF_0002.1", "1", "Species alpha", "diversity", "True", 90),
                inventory_row("GCF_0003.1", "2", "Species beta", "canonical", "False", 100),
                inventory_row("GCF_0004.1", "3", "Species gamma", "canonical", "True", 100),
                inventory_row("GCF_0005.1", "4", "Bacillus cereus", "canonical", "False", 100),
            ]
            write_tsv(inventory, inventory_rows)
            resolution = root / "resolution.tsv"
            write_tsv(
                resolution,
                [
                    resolution_row("GCF_0001.1", "", "missing"),
                    resolution_row("GCF_0002.1", "GCF_0002.1", "exact"),
                    resolution_row(
                        "GCF_0003.1", "GCF_0003.2", "current_version_substituted"
                    ),
                    resolution_row("GCF_0004.1", "", "missing"),
                    resolution_row("GCF_0005.1", "GCF_0005.1", "exact"),
                ],
            )
            reports = root / "reports.jsonl"
            reports.write_text(
                json.dumps(report("GCF_0002.1", "1", "Species alpha"))
                + "\n"
                + json.dumps(report("GCF_0003.2", "22", "Species beta revised"))
                + "\n"
                + json.dumps(report("GCF_0005.1", "4", "Bacillus cereus group sp. N1"))
                + "\n"
            )
            pre_missing = root / "pre_missing.tsv"
            write_tsv(
                pre_missing,
                [
                    {
                        "taxid": "99",
                        "species_name": "Species absent",
                        "exclusion_reasons": "no_matching_gtdb_metadata",
                    }
                ],
            )

            summary = reconcile_reference_inventory(
                inventory,
                resolution,
                reports,
                root / "output",
                pre_inventory_missing_tsv=pre_missing,
            )

            self.assertEqual(summary["current_reference_count"], 1)
            self.assertEqual(summary["missing_species_count"], 3)
            self.assertEqual(summary["missing_controlled_species_count"], 1)
            self.assertEqual(summary["total_controlled_species_gap_count"], 2)
            self.assertEqual(summary["taxonomy_review_species_count"], 1)
            self.assertEqual(summary["excluded_current_metadata_count"], 2)
            self.assertEqual(
                summary["excluded_current_metadata_reason_counts"],
                {
                    "taxonomy_assignment_conflict": 1,
                    "unnamed_current_species": 1,
                },
            )
            self.assertEqual(summary["taxonomy_review_reference_count"], 1)
            current = list(
                csv.DictReader(
                    (root / "output/reference_inventory_current.tsv").open(),
                    delimiter="\t",
                )
            )
            by_accession = {row["accession"]: row for row in current}
            self.assertEqual(by_accession["GCF_0002.1"]["reference_role"], "canonical")
            self.assertIn("Promoted", by_accession["GCF_0002.1"]["reason"])
            self.assertNotIn("GCF_0003.2", by_accession)
            gaps = (root / "output/controlled_species_gaps.tsv").read_text()
            self.assertIn("gtdb_inventory\t99\tSpecies absent", gaps)
            self.assertIn("ncbi_current_resolution\t3\tSpecies gamma", gaps)
            review_rows = list(
                csv.DictReader(
                    (root / "output/taxonomy_review.tsv").open(), delimiter="\t"
                )
            )
            self.assertEqual([row["accession"] for row in review_rows], ["GCF_0003.2"])

    def test_official_merge_and_direct_subspecies_do_not_create_false_conflicts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inventory = root / "inventory.tsv"
            merged_alias = inventory_row(
                "GCF_0010.1", "20", "Species current", "canonical", "True", 100
            )
            merged_alias.update(source_species_taxid="10", source_taxid="10")
            subspecies = inventory_row(
                "GCF_0040.1", "40", "Species parent", "canonical", "True", 100
            )
            subspecies.update(source_species_taxid="30", source_taxid="40")
            write_tsv(inventory, [merged_alias, subspecies])
            resolution = root / "resolution.tsv"
            write_tsv(
                resolution,
                [
                    resolution_row("GCF_0010.1", "GCF_0010.1", "exact"),
                    resolution_row("GCF_0040.1", "GCF_0040.1", "exact"),
                ],
            )
            reports = root / "reports.jsonl"
            reports.write_text(
                json.dumps(report("GCF_0010.1", "20", "Species current"))
                + "\n"
                + json.dumps(
                    report(
                        "GCF_0040.1",
                        "30",
                        "Species parent",
                        organism_taxid="40",
                    )
                )
                + "\n"
            )
            taxonomy = root / "taxonomy"
            taxonomy.mkdir()
            (taxonomy / "merged.dmp").write_text("10\t|\t20\t|\n")

            summary = reconcile_reference_inventory(
                inventory,
                resolution,
                reports,
                root / "output",
                ncbi_taxonomy_dir=taxonomy,
            )

            self.assertEqual(summary["current_reference_count"], 2)
            self.assertEqual(summary["taxonomy_review_reference_count"], 0)
            self.assertEqual(summary["excluded_current_metadata_count"], 0)


def inventory_row(
    accession: str,
    taxid: str,
    species: str,
    role: str,
    controlled: str,
    quality: int,
) -> dict[str, str]:
    return {
        "accession": accession,
        "species_taxid": taxid,
        "ncbi_species_name": species,
        "gtdb_species_name": species,
        "reference_role": role,
        "controlled_target": controlled,
        "quality_score": str(quality),
        "reason": "Selected for test.",
    }


def resolution_row(requested: str, returned: str, status: str) -> dict[str, str]:
    return {
        "requested_accession": requested,
        "returned_accession": returned,
        "status": status,
    }


def report(
    accession: str,
    taxid: str,
    species: str,
    *,
    organism_taxid: str | None = None,
) -> dict[str, object]:
    return {
        "accession": accession,
        "assembly_info": {"assembly_status": "current"},
        "organism": {"tax_id": organism_taxid or taxid, "organism_name": species},
        "checkm_info": {"checkm_species_tax_id": taxid},
        "average_nucleotide_identity": {
            "submitted_species": species,
            "taxonomy_check_status": "OK",
        },
    }


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fields = sorted({field for row in rows for field in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
