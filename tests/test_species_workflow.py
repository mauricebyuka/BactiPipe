from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from bactipipe.species_identification import AniHit
from bactipipe.speciesfinder_database import SpeciesFinderDatabase
from bactipipe.speciesfinder_evidence import SpeciesFinderHit
from bactipipe.species_workflow import identify_assembly


class SpeciesWorkflowTests(unittest.TestCase):
    def test_speciesfinder_no_candidates_fails_closed_without_crashing(self):
        self.check_speciesfinder_input("Illumina", 0)

    def test_speciesfinder_uses_assembly_even_when_nanopore_reads_available(self):
        self.check_speciesfinder_input("Nanopore", 1)

    def test_speciesfinder_uses_assembly_even_when_illumina_pairs_available(self):
        self.check_speciesfinder_input("Illumina", 2)

    def test_failed_nomination_preserves_previous_results_and_original(self):
        self.check_speciesfinder_input("Nanopore", 1, fail=True)

    def check_speciesfinder_input(self, platform, read_count, fail=False):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            assembly = root / "S1.fasta"
            assembly.write_text(">S1\n" + "ACGT" * 10001 + "\n")
            original = assembly.read_bytes()
            reads = tuple(root / f"reads{index}.fastq" for index in range(read_count))
            for read in reads:
                read.write_text("@read\nACGT\n+\nIIII\n")
            database = root / "database"
            database.mkdir()
            output = root / "results" / "S1"
            identity = {
                "version": "speciesfinder-candidate-v1",
                "status": "candidate",
                "database_kind": "speciesfinder_kma",
            }
            (database / "database_identity.json").write_text(json.dumps(identity))
            prepared = SpeciesFinderDatabase(
                root=database,
                database_prefix=database / "bacteria",
                taxonomy_file=database / "bacteria.tax",
                template_lookup=database / "lookup.sqlite3",
                identity=identity,
            )

            def fake_nominate(*, inputs, **kwargs):
                self.assertNotEqual(inputs, (assembly,))
                query = Path(inputs[0])
                self.assertEqual(query.name, "query.fasta")
                self.assertEqual([len(line) for line in query.read_text().splitlines()
                                  if not line.startswith(">")], [20000, 20000, 4])
                self.assertTrue((query.parent / "query.json").is_file())
                if fail:
                    raise RuntimeError("test nomination failure")
                return ()

            def fake_extract(*, output_directory, **_kwargs):
                candidate_root = Path(output_directory)
                candidate_root.mkdir(parents=True, exist_ok=True)
                (candidate_root / "reference_manifest.tsv").write_text(
                    "accession\tspecies_taxid\tspecies_name\tfasta_path\tsource\t"
                    "reference_role\tassembly_level\tsequence_origin\ttype_material\t"
                    "canonical_for_tree\tcallable_for_identification\tinclude_in_kraken\t"
                    "reporting_group\tncbi_taxonomy_name\tgtdb_taxonomy\tsha256\n"
                )
                return ()

            with (
                patch(
                    "bactipipe.species_workflow.validate_speciesfinder_database",
                    return_value=prepared,
                ),
                patch(
                    "bactipipe.species_workflow.run_speciesfinder",
                    side_effect=fake_nominate,
                ) as nominate,
                patch(
                    "bactipipe.species_workflow.extract_candidate_references",
                    side_effect=fake_extract,
                ),
                patch(
                    "bactipipe.species_workflow.command_version",
                    return_value="v1.6.13",
                ),
                patch(
                    "bactipipe.species_workflow.python_package_version",
                    return_value="v4.0.4",
                ),
            ):
                if fail:
                    output.mkdir(parents=True)
                    (output / "previous.txt").write_text("prior complete result")
                    with self.assertRaisesRegex(RuntimeError, "test nomination failure"):
                        identify_assembly(sample_id="S1", assembly=assembly,
                                          database_root=database, output_directory=output,
                                          platform=platform, reads=reads, allow_candidate_database=True)
                    self.assertEqual((output / "previous.txt").read_text(), "prior complete result")
                    self.assertEqual(assembly.read_bytes(), original)
                    self.assertEqual(list(output.parent.glob(".*.building-*")), [])
                    return
                output.mkdir(parents=True)
                (output / "previous_intact_nomination.json").write_text('{"schema_version": 3}')
                result = identify_assembly(
                    sample_id="S1",
                    assembly=assembly,
                    database_root=database,
                    output_directory=output,
                    platform=platform,
                    reads=reads,
                    allow_candidate_database=True,
                    require_mixed_signal_assessment=True,
                )

            self.assertEqual(assembly.read_bytes(), original)
            self.assertFalse((output / "previous_intact_nomination.json").exists())
            self.assertFalse((output / "speciesfinder_query/query.fasta").exists())
            self.assertTrue((output / "speciesfinder_query/fragments.tsv").exists())
            provenance = json.loads((output / "analysis_provenance.json").read_text())
            self.assertEqual(provenance["parameters"]["speciesfinder_input_kind"], "fragmented_assembly_fasta")
            self.assertEqual(provenance["schema_version"], 4)
            self.assertEqual(provenance["parameters"]["speciesfinder_query"]["limit_bp"], 20000)
            self.assertEqual(len(provenance["reads"]), read_count)
            self.assertEqual(result.status, "mixed_signal_review")
            self.assertEqual(result.reported_organism, "No confident identification")
            self.assertEqual(
                (output / "taxonomic_signals.tsv").read_text(),
                "species_taxid\tspecies_name\tread_count\tread_percent\t"
                "distinct_minimizers\treference_breadth\tcompetitive_mapped_reads\t"
                "competitive_mapped_percent\tnomination_score\t"
                "nomination_percent\n",
            )

    def test_speciesfinder_nominates_candidates_but_ani_makes_the_call(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            assembly = root / "S1.fasta"
            assembly.write_text(">S1\nACGTACGT\n")
            database = root / "database"
            database.mkdir()
            output = root / "results" / "S1"
            identity = {
                "version": "speciesfinder-candidate-v1",
                "status": "candidate",
                "database_kind": "speciesfinder_kma",
            }
            (database / "database_identity.json").write_text(json.dumps(identity))
            prepared = SpeciesFinderDatabase(
                root=database,
                database_prefix=database / "bacteria",
                taxonomy_file=database / "bacteria.tax",
                template_lookup=database / "lookup.sqlite3",
                identity=identity,
            )
            nomination = SpeciesFinderHit(
                assembly_accession="GCF_1",
                template_accession="NC_1.1",
                description="Species alpha",
                organism_taxid="1",
                taxonomy="Bacteria; Species alpha",
                species_taxid="1",
                species_name="Species alpha",
                query_count=0,
                score=100,
                query_identity=99,
                query_coverage=95,
                template_coverage=95,
                depth=1,
            )
            ani = AniHit(
                accession="NC_1.1__kma1",
                species_taxid="1",
                species_name="Species alpha",
                reference_path="NC_1.1__kma1.fna",
                ani=99.0,
                aligned_fraction_reference=95.0,
                aligned_fraction_query=95.0,
            )

            def fake_extract(*, output_directory, **_kwargs):
                candidate_root = Path(output_directory)
                candidate_root.mkdir(parents=True, exist_ok=True)
                reference = candidate_root / "reference.fna"
                reference.write_text(">NC_1.1\nACGT\n")
                (candidate_root / "reference_manifest.tsv").write_text(
                    "accession\tspecies_taxid\tspecies_name\tfasta_path\tsource\t"
                    "reference_role\tassembly_level\tsequence_origin\ttype_material\t"
                    "canonical_for_tree\tcallable_for_identification\tinclude_in_kraken\t"
                    "reporting_group\tncbi_taxonomy_name\tgtdb_taxonomy\tsha256\n"
                    "NC_1.1__kma1\t1\tSpecies alpha\treference.fna\tSpeciesFinder\t"
                    "representative\tcomplete\tdatabase_reference\tno\tyes\tyes\tno\t"
                    "\tSpecies alpha\t\tchecksum\n"
                )
                return (MagicMock(fasta_path=reference),)

            def fake_dist(*, output_path, **_kwargs):
                self.assertEqual(_kwargs["assembly"], assembly)
                Path(output_path).write_text("mock\n")

            with (
                patch(
                    "bactipipe.species_workflow.validate_speciesfinder_database",
                    return_value=prepared,
                ),
                patch(
                    "bactipipe.species_workflow.run_speciesfinder",
                    return_value=(nomination,),
                ),
                patch(
                    "bactipipe.species_workflow.extract_candidate_references",
                    side_effect=fake_extract,
                ),
                patch(
                    "bactipipe.species_workflow.run_skani_candidate_dist",
                    side_effect=fake_dist,
                ),
                patch(
                    "bactipipe.species_workflow.parse_skani_output",
                    return_value=(ani,),
                ),
                patch(
                    "bactipipe.species_workflow.command_version",
                    return_value="v1.6.13",
                ),
                patch(
                    "bactipipe.species_workflow.python_package_version",
                    return_value="v4.0.4",
                ),
            ):
                result = identify_assembly(
                    sample_id="S1",
                    assembly=assembly,
                    database_root=database,
                    output_directory=output,
                    platform="Illumina",
                    allow_candidate_database=True,
                )

            self.assertEqual(result.status, "confirmed_species")
            self.assertEqual(result.reported_organism, "Species alpha")
            provenance = json.loads((output / "analysis_provenance.json").read_text())
            self.assertEqual(
                provenance["parameters"]["candidate_provider"], "speciesfinder"
            )
            self.assertEqual(provenance["tools"]["speciesfinder"], "v4.0.4")

    def test_required_missing_read_evidence_fails_closed_and_records_provenance(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            assembly = root / "S1.fasta"
            assembly.write_text(">S1\nACGTACGT\n")
            database = root / "database"
            database.mkdir()
            output = root / "results" / "S1"
            hit = AniHit(
                accession="GCF_1",
                species_taxid="1",
                species_name="Species alpha",
                reference_path="GCF_1.fna",
                ani=99.0,
                aligned_fraction_reference=95.0,
                aligned_fraction_query=95.0,
            )

            def fake_search(*, output_path, **_kwargs):
                Path(output_path).write_text("mock\n")

            with (
                patch(
                    "bactipipe.species_workflow.validate_database",
                    return_value={"version": "candidate-v1", "status": "candidate"},
                ),
                patch(
                    "bactipipe.species_workflow.ReferenceCatalog.from_tsv",
                    return_value=MagicMock(),
                ),
                patch(
                    "bactipipe.species_workflow.run_skani_search",
                    side_effect=fake_search,
                ),
                patch(
                    "bactipipe.species_workflow.parse_skani_output",
                    return_value=(hit,),
                ),
                patch(
                    "bactipipe.species_workflow.command_version",
                    return_value="v0.3.0",
                ),
            ):
                result = identify_assembly(
                    sample_id="S1",
                    assembly=assembly,
                    database_root=database,
                    output_directory=output,
                    platform="Illumina",
                    require_mixed_signal_assessment=True,
                    allow_candidate_database=True,
                )

            self.assertEqual(result.status, "mixed_signal_review")
            self.assertIsNone(result.match_percent)
            payload = json.loads((output / "analysis_provenance.json").read_text())
            self.assertEqual(payload["database"]["version"], "candidate-v1")
            self.assertEqual(payload["tools"]["skani"], "v0.3.0")
            self.assertTrue(
                payload["parameters"]["mixed_signal_assessment_required"]
            )
            self.assertEqual(payload["assembly"]["filename"], "S1.fasta")


if __name__ == "__main__":
    unittest.main()
