from __future__ import annotations

import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path

from bactipipe.reference_curation import (
    CurationPolicy,
    plan_gtdb_reference_inventory,
    plan_references,
    read_merged_taxids,
    write_inventory_plan,
)


def report(
    accession: str,
    *,
    species_taxid: int = 75985,
    species: str = "Mannheimia haemolytica",
    completeness: float = 99.0,
    contamination: float = 1.0,
    type_material: bool = False,
    refseq_category: str = "",
    model: str = "Generic",
) -> dict:
    value = {
        "accession": accession,
        "source_database": "SOURCE_DATABASE_REFSEQ",
        "organism": {"organism_name": species, "tax_id": species_taxid},
        "assembly_info": {
            "assembly_status": "current",
            "assembly_level": "Complete Genome",
            "refseq_category": refseq_category,
            "biosample": {"models": [model]},
        },
        "assembly_stats": {
            "number_of_contigs": 1,
            "contig_n50": 2_500_000,
            "total_sequence_length": 2_500_000,
        },
        "checkm_info": {
            "checkm_species_tax_id": species_taxid,
            "completeness": completeness,
            "contamination": contamination,
        },
        "average_nucleotide_identity": {
            "submitted_species": species,
            "submitted_organism": species,
            "taxonomy_check_status": "OK",
        },
    }
    if type_material:
        value["type_material"] = {"type_label": "TYPE_MATERIAL"}
    return value


class ReferenceCurationTests(unittest.TestCase):
    def test_invalid_policy_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "minimum_completeness"):
            CurationPolicy(minimum_completeness=101)

    def test_type_material_wins_over_non_type_reference(self):
        planned = plan_references(
            [
                report("GCF_TYPE", type_material=True),
                report("GCF_REF", refseq_category="reference genome"),
            ]
        )
        decisions = {row.accession: row.decision for row in planned}
        self.assertEqual(decisions["GCF_TYPE"], "selected_canonical")
        self.assertEqual(decisions["GCF_REF"], "eligible_not_selected")

    def test_failed_quality_and_mag_are_audited(self):
        planned = plan_references(
            [
                report("GCF_LOW", completeness=70.0),
                report("GCF_MAG", model="Metagenome-assembled Genome"),
            ]
        )
        reasons = {row.accession: row.reason for row in planned}
        self.assertIn("low_or_missing_completeness", reasons["GCF_LOW"])
        self.assertIn("non_isolate_source", reasons["GCF_MAG"])
        self.assertTrue(all(row.decision == "excluded" for row in planned))

    def test_species_are_selected_independently(self):
        planned = plan_references(
            [
                report("GCF_1", species_taxid=1, species="Species alpha"),
                report("GCF_2", species_taxid=2, species="Species beta"),
            ]
        )
        self.assertEqual(
            sum(row.decision == "selected_canonical" for row in planned), 2
        )

    def test_gtdb_inventory_is_bounded_to_targets_and_same_genus_representatives(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            catalog = root / "controlled.txt"
            catalog.write_text(
                "75985\t2500000\tNZ_TEST\tMannheimia haemolytica\n"
                "789\t2400000\tNZ_MAG\tMannheimia environmentalensis\n"
            )
            metadata = root / "metadata.tsv.gz"
            taxonomy = taxonomy_dir(root)
            rows = [
                gtdb_row(
                    "RS_GCF_0001.1",
                    taxid="75985",
                    species="Mannheimia haemolytica",
                    representative="f",
                    type_material="assembly from type material",
                ),
                gtdb_row(
                    "RS_GCF_0002.1",
                    taxid="123",
                    species="Mannheimia glucosida",
                    representative="t",
                ),
                gtdb_row(
                    "RS_GCF_0003.1",
                    taxid="456",
                    species="Pasteurella multocida",
                    representative="t",
                ),
                gtdb_row(
                    "GB_GCA_0004.1",
                    taxid="789",
                    species="Mannheimia environmentalensis",
                    representative="t",
                    genome_category="derived from metagenome",
                ),
                gtdb_row(
                    "GB_GCA_0005.1",
                    taxid="999",
                    species="Mannheimia sp000005",
                    representative="t",
                    ncbi_species="",
                ),
                gtdb_row(
                    "RS_GCF_0006.1",
                    taxid="998",
                    species="Mannheimia haemolytica",
                    representative="t",
                ),
            ]
            with gzip.open(metadata, "wt", newline="") as handle:
                writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)

            references, summary = plan_gtdb_reference_inventory(
                metadata, catalog, ncbi_taxonomy_dir=taxonomy
            )

            self.assertEqual(
                {record.accession for record in references},
                {"GCF_0001.1", "GCF_0002.1", "GCF_0006.1"},
            )
            roles = {record.accession: record.reference_role for record in references}
            self.assertEqual(roles["GCF_0001.1"], "canonical")
            self.assertEqual(roles["GCF_0002.1"], "canonical")
            by_accession = {record.accession: record for record in references}
            self.assertFalse(by_accession["GCF_0006.1"].controlled_target)
            self.assertEqual(summary["counts"]["excluded:non_isolate_source"], 1)
            self.assertEqual(
                summary["counts"]["excluded:unresolved_ncbi_species_name"], 1
            )
            self.assertEqual(summary["counts"]["controlled_species_missing"], 1)
            self.assertEqual(
                summary["missing_controlled_species"][0]["exclusion_reasons"],
                "non_isolate_source=1",
            )

    def test_inventory_outputs_are_atomic_and_machine_readable(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            catalog = root / "controlled.txt"
            catalog.write_text("75985\t2500000\tNZ_TEST\tMannheimia haemolytica\n")
            metadata = root / "metadata.tsv"
            taxonomy = taxonomy_dir(root)
            row = gtdb_row(
                "RS_GCF_0001.1",
                taxid="75985",
                species="Mannheimia haemolytica",
                representative="t",
            )
            with metadata.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(row))
                writer.writeheader()
                writer.writerow(row)
            references, summary = plan_gtdb_reference_inventory(
                metadata, catalog, ncbi_taxonomy_dir=taxonomy
            )
            write_inventory_plan(
                references,
                output_tsv=root / "inventory.tsv",
                accessions_path=root / "accessions.txt",
                summary=summary,
                summary_path=root / "summary.json",
                missing_controlled_path=root / "missing.tsv",
                controlled_aliases_path=root / "aliases.tsv",
            )
            self.assertEqual((root / "accessions.txt").read_text(), "GCF_0001.1\n")
            self.assertEqual(json.loads((root / "summary.json").read_text()), summary)
            self.assertEqual(
                (root / "missing.tsv").read_text().splitlines()[0],
                "taxid\texpected_genome_size\tlegacy_sequence_accession\tspecies_name\tmetadata_rows\texclusion_reasons",
            )

    def test_merged_taxids_and_direct_subspecies_targets_are_audited(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            catalog = root / "controlled.txt"
            catalog.write_text(
                "10\t2000000\tNZ_OLD\tSpecies old name\n"
                "20\t2000000\tNZ_CURRENT\tSpecies current name\n"
                "40\t2000000\tNZ_SUB\tStreptococcus equi subsp zooepidemicus\n"
            )
            taxonomy = taxonomy_dir(root, "10\t|\t20\t|\n")
            metadata = root / "metadata.tsv"
            rows = [
                gtdb_row(
                    "RS_GCF_0010.1",
                    taxid="20",
                    species="Species current name",
                    representative="t",
                ),
                gtdb_row(
                    "RS_GCF_0040.1",
                    taxid="30",
                    source_taxid="40",
                    species="Streptococcus equi",
                    representative="t",
                ),
            ]
            with metadata.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)

            references, summary = plan_gtdb_reference_inventory(
                metadata, catalog, ncbi_taxonomy_dir=taxonomy
            )

            by_accession = {record.accession: record for record in references}
            merged = by_accession["GCF_0010.1"]
            self.assertEqual(merged.species_taxid, "20")
            self.assertEqual(merged.controlled_catalog_taxid, "10;20")
            self.assertEqual(merged.controlled_match, "merged_ncbi_taxid")
            subspecies = by_accession["GCF_0040.1"]
            self.assertEqual(subspecies.species_taxid, "40")
            self.assertEqual(subspecies.source_species_taxid, "30")
            self.assertEqual(subspecies.source_taxid, "40")
            self.assertEqual(subspecies.controlled_match, "ncbi_taxid")
            self.assertEqual(summary["counts"]["controlled_species_missing"], 0)
            self.assertEqual(summary["counts"]["controlled_taxid_alias_group_count"], 1)
            self.assertEqual(
                summary["controlled_taxid_aliases"][0]["catalog_taxids"], "10;20"
            )

    def test_merged_taxid_reader_rejects_cycles_when_resolved(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            taxonomy = taxonomy_dir(root, "10\t|\t20\t|\n20\t|\t10\t|\n")
            merged = read_merged_taxids(taxonomy)
            catalog = root / "controlled.txt"
            catalog.write_text("10\t2000000\tNZ_TEST\tSpecies alpha\n")
            metadata = root / "metadata.tsv"
            row = gtdb_row(
                "RS_GCF_0010.1", taxid="10", species="Species alpha", representative="t"
            )
            with metadata.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(row))
                writer.writeheader()
                writer.writerow(row)
            with self.assertRaisesRegex(ValueError, "Cycle detected"):
                plan_gtdb_reference_inventory(
                    metadata, catalog, ncbi_taxonomy_dir=taxonomy
                )


def gtdb_row(
    accession: str,
    *,
    taxid: str,
    source_taxid: str | None = None,
    species: str,
    representative: str,
    type_material: str = "none",
    genome_category: str = "none",
    ncbi_species: str | None = None,
) -> dict[str, str]:
    genus = species.split()[0]
    ncbi_species = species if ncbi_species is None else ncbi_species
    return {
        "accession": accession,
        "gtdb_representative": representative,
        "gtdb_taxonomy": f"d__Bacteria;p__P;c__C;o__O;f__F;g__{genus};s__{species}",
        "ncbi_species_taxid": taxid,
        "ncbi_taxid": source_taxid or taxid,
        "ncbi_taxonomy": f"d__Bacteria;p__P;c__C;o__O;f__F;g__{genus};s__{ncbi_species}",
        "ncbi_organism_name": f"{species} isolate test",
        "ncbi_genome_category": genome_category,
        "ncbi_genome_representation": "full",
        "ncbi_assembly_level": "Complete Genome",
        "ncbi_refseq_category": "representative genome",
        "ncbi_type_material_designation": type_material,
        "checkm2_completeness": "99.5",
        "checkm2_contamination": "0.2",
        "checkm_completeness": "99.0",
        "checkm_contamination": "0.4",
        "ncbi_contig_count": "none",
        "ncbi_contig_n50": "none",
        "ncbi_total_length": "none",
        "contig_count": "1",
        "n50_contigs": "2500000",
        "genome_size": "2500000",
    }


def taxonomy_dir(root: Path, merged: str = "") -> Path:
    taxonomy = root / "taxonomy"
    taxonomy.mkdir(exist_ok=True)
    (taxonomy / "merged.dmp").write_text(merged)
    return taxonomy


if __name__ == "__main__":
    unittest.main()
