from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from bactipipe.reference_supplements import augment_inventory_with_gap_candidates


class ReferenceSupplementTests(unittest.TestCase):
    def test_only_selected_exact_controlled_gap_is_added(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inventory = root / "inventory.tsv"
            write_tsv(inventory, [{"accession": "GCF_0001.1", "species_taxid": "1", "ncbi_species_name": "Species alpha"}])
            candidates = root / "candidates.tsv"
            write_tsv(
                candidates,
                [
                    {
                        "accession": "GCF_0002.1",
                        "species_taxid": "2",
                        "species_name": "Species beta",
                        "submitted_organism": "Species beta strain T",
                        "source": "REFSEQ",
                        "assembly_level": "Complete Genome",
                        "refseq_category": "reference genome",
                        "type_material": "True",
                        "sequence_origin": "isolate",
                        "completeness": "99",
                        "contamination": "0",
                        "contigs": "1",
                        "n50": "2000000",
                        "total_length": "2000000",
                        "taxonomy_check_status": "OK",
                        "quality_score": "10",
                        "decision": "selected_canonical",
                    }
                ],
            )
            catalog = root / "controlled.txt"
            catalog.write_text("2\t2000000\tNZ_TEST\tSpecies beta\n")
            gaps = root / "gaps.tsv"
            write_tsv(
                gaps,
                [{"stage": "test", "species_taxid": "2", "species_name": "Species beta", "reason": "missing"}],
            )
            summary = augment_inventory_with_gap_candidates(
                inventory, candidates, catalog, gaps, root / "out"
            )
            self.assertEqual(summary["supplement_reference_count"], 1)
            self.assertEqual(summary["controlled_gap_count_after_recovery"], 0)
            rows = list(csv.DictReader((root / "out/reference_inventory_augmented.tsv").open(), delimiter="\t"))
            supplement = next(row for row in rows if row["accession"] == "GCF_0002.1")
            self.assertEqual(supplement["controlled_catalog_name"], "Species beta")
            self.assertEqual(supplement["reference_role"], "canonical")


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fields = sorted({field for row in rows for field in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
