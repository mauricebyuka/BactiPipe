from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from bactipipe.speciesfinder_database import (
    audit_speciesfinder_taxonomy,
    extract_candidate_references,
    prepare_speciesfinder_database,
    validate_speciesfinder_database,
)
from bactipipe.speciesfinder_evidence import SpeciesFinderHit


class SpeciesFinderDatabaseTests(unittest.TestCase):
    def test_taxonomy_audit_excludes_mag_and_unresolved_records(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            taxonomy = root / "bacteria.tax"
            taxonomy.write_text(
                "R1 Species alpha complete genome\tGCA_1\t1\tBacteria\t1\tSpecies alpha\n"
                "R2 MAG: Species beta isolate MAG_1\tGCA_2\t2\tBacteria\t2\tSpecies beta\n"
                "R3 Species gamma complete genome\tNA\tNA\tNA\tNA\tSpecies gamma\n"
                "R4 unknown bacterial sequence\tNA\tNA\tNA\tNA\tNA\n"
                "R5 Enterobacter isolate\tGCA_5\t5\tBacteria\t5\tEnterobacter sp.\n"
            )
            audit = root / "audit.tsv"
            counts = audit_speciesfinder_taxonomy(taxonomy, audit)
            self.assertEqual(counts["eligible_nomination_template"], 1)
            self.assertEqual(counts["excluded_non_isolate_source"], 1)
            self.assertEqual(counts["excluded_unresolved_species_taxid"], 1)
            self.assertEqual(counts["excluded_unresolved_species_name"], 1)
            self.assertEqual(counts["excluded_unresolved_species_label"], 1)
            self.assertEqual(len(audit.read_text().splitlines()), 6)

    def test_prepare_is_immutable_and_candidate_by_default(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prefix, taxonomy, checksum = self._source(root)
            with patch(
                "bactipipe.speciesfinder_database._tool_version",
                return_value="KMA-1.6.13",
            ):
                result = prepare_speciesfinder_database(
                    database_prefix=prefix,
                    taxonomy_file=taxonomy,
                    archive_sha256_file=checksum,
                    output_root=root / "prepared",
                    version="test-v1",
                )
            database = validate_speciesfinder_database(
                result, allow_candidate=True
            )
            self.assertEqual(database.identity["database_kind"], "speciesfinder_kma")
            self.assertEqual(database.identity["template_count"], 2)
            with self.assertRaisesRegex(ValueError, "not an active production"):
                validate_speciesfinder_database(result)

    def test_extracts_only_nominated_template_and_writes_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prefix, taxonomy, checksum = self._source(root)
            with patch(
                "bactipipe.speciesfinder_database._tool_version",
                return_value="KMA-1.6.13",
            ):
                prepared = prepare_speciesfinder_database(
                    database_prefix=prefix,
                    taxonomy_file=taxonomy,
                    archive_sha256_file=checksum,
                    output_root=root / "prepared",
                    version="test-v1",
                )
            database = validate_speciesfinder_database(
                prepared, allow_candidate=True
            )
            hit = SpeciesFinderHit(
                assembly_accession="GCF_1",
                template_accession="NC_1.1",
                description="Species alpha chromosome",
                organism_taxid="1",
                taxonomy="Bacteria; Species alpha",
                species_taxid="1",
                species_name="Species alpha",
                query_count=100,
                score=200,
                query_identity=99,
                query_coverage=98,
                template_coverage=97,
                depth=10,
            )

            def fake_run(_command, **kwargs):
                kwargs["stdout"].write(">NC_1.1 Species alpha chromosome\nACGTACGT\n")
                return SimpleNamespace(returncode=0, stderr="")

            with patch(
                "bactipipe.speciesfinder_database.subprocess.run",
                side_effect=fake_run,
            ):
                records = extract_candidate_references(
                    database=database,
                    hits=[hit],
                    output_directory=root / "sample_candidates",
                )
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].template_id, 1)
            self.assertTrue(records[0].fasta_path.is_file())
            manifest = root / "sample_candidates/reference_manifest.tsv"
            self.assertIn("NC_1.1__kma1", manifest.read_text())

    @staticmethod
    def _source(root: Path) -> tuple[Path, Path, Path]:
        prefix = root / "bacteria"
        for suffix in (".comp.b", ".seq.b", ".length.b"):
            Path(f"{prefix}{suffix}").write_bytes(b"index")
        Path(f"{prefix}.name").write_text(
            "NC_1.1 Species alpha chromosome\nNC_2.1 Species beta chromosome\n"
        )
        taxonomy = root / "bacteria.tax"
        taxonomy.write_text("NC_1.1 Species alpha\tGCF_1\t1\tBacteria\t1\tSpecies alpha\n")
        checksum = root / "archive.sha256"
        checksum.write_text(hashlib.sha256(b"archive").hexdigest() + "  archive\n")
        return prefix, taxonomy, checksum


if __name__ == "__main__":
    unittest.main()
