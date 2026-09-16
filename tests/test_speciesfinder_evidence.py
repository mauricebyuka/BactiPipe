from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from bactipipe.speciesfinder_evidence import (
    parse_speciesfinder_results,
    run_speciesfinder,
    select_candidate_hits,
    speciesfinder_hits_to_signals,
    is_reportable_species_name,
)


class SpeciesFinderEvidenceTests(unittest.TestCase):
    def test_parser_preserves_metrics_and_controlled_subspecies(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = Path(temporary) / "results.txt"
            result.write_text(
                "# Assembly\tNum\tScore\tTemplate_length\tQuery_Identity\t"
                "Query_Coverage\tTemplate_Coverage\tDepth\tAccession Number\t"
                "Description\tTAXID\tTaxonomy\tTAXID\tSpecies\n"
                "GCF_1\t800\t900\t2100000\t99.3\t98.0\t97.0\t20\t"
                "NZ_CP1.1\tStreptococcus equi subsp. zooepidemicus genome\t"
                "40041\tBacteria; Streptococcus equi subsp. zooepidemicus\t"
                "1336\tStreptococcus equi\n"
            )
            hits = parse_speciesfinder_results(result)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].species_taxid, "40041")
        self.assertEqual(
            hits[0].species_name,
            "Streptococcus equi subsp. zooepidemicus",
        )
        self.assertEqual(hits[0].query_count, 800)
        self.assertEqual(hits[0].query_identity, 99.3)
        self.assertEqual(hits[0].organism_taxid, "40041")

    def test_candidate_selection_is_bounded_per_species(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = Path(temporary) / "results.txt"
            result.write_text(
                "# Assembly\tScore\tQuery_Coverage\tTemplate_Coverage\t"
                "Accession Number\tDescription\tTAXID\tTaxonomy\t"
                "TAXID Species\tSpecies\n"
                "A\t100\t99\t99\tR1\ta\t1\tBacteria\t1\tSpecies alpha\n"
                "A\t90\t98\t98\tR2\ta\t1\tBacteria\t1\tSpecies alpha\n"
                "B\t80\t97\t97\tR3\tb\t2\tBacteria\t2\tSpecies beta\n"
            )
            hits = parse_speciesfinder_results(result)
        selected = select_candidate_hits(
            hits, maximum_species=2, references_per_species=1
        )
        self.assertEqual([hit.template_accession for hit in selected], ["R1", "R3"])

    def test_query_coverage_above_one_hundred_is_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = Path(temporary) / "results.txt"
            result.write_text(
                "# Assembly\tScore\tQuery_Identity\tQuery_Coverage\t"
                "Template_Coverage\tAccession Number\tDescription\tTAXID\t"
                "Taxonomy\tTAXID Species\tSpecies\n"
                "A\t100\t99.9\t1581.92\t78.69\tR1\ta\t1\tBacteria\t"
                "1\tSpecies alpha\n"
            )
            hits = parse_speciesfinder_results(result)
        self.assertEqual(hits[0].query_coverage, 1581.92)

    def test_template_coverage_above_one_hundred_is_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = Path(temporary) / "results.txt"
            result.write_text(
                "# Assembly\tScore\tQuery_Identity\tQuery_Coverage\t"
                "Template_Coverage\tAccession Number\tDescription\tTAXID\t"
                "Taxonomy\tTAXID Species\tSpecies\n"
                "A\t100\t99.9\t99.99\t100.01\tR1\ta\t1\tBacteria\t"
                "1\tSpecies alpha\n"
            )
            hits = parse_speciesfinder_results(result)
        self.assertEqual(hits[0].query_coverage, 99.99)
        self.assertEqual(hits[0].template_coverage, 100.01)

    def test_score_ranks_candidates_before_length_ratio(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = Path(temporary) / "results.txt"
            result.write_text(
                "# Assembly\tScore\tQuery_Coverage\tTemplate_Coverage\t"
                "Accession Number\tDescription\tTAXID\tTaxonomy\t"
                "TAXID Species\tSpecies\n"
                "A\t100\t1500\t5\tR1\ta\t1\tBacteria\t1\tSpecies alpha\n"
                "B\t1000\t100\t90\tR2\tb\t2\tBacteria\t2\tSpecies beta\n"
            )
            hits = parse_speciesfinder_results(result)
        self.assertEqual(hits[0].template_accession, "R2")

    def test_metagenome_assembled_template_cannot_be_a_candidate(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = Path(temporary) / "results.txt"
            result.write_text(
                "# Assembly\tScore\tQuery_Coverage\tTemplate_Coverage\t"
                "Accession Number\tDescription\tTAXID\tTaxonomy\tTAXID\tSpecies\n"
                "A\t1000\t100\t99\tR1\tMAG: Species alpha isolate MAG_1\t"
                "1\tBacteria\t1\tSpecies alpha\n"
                "B\t900\t100\t98\tR2\tSpecies beta isolate culture 1\t"
                "2\tBacteria\t2\tSpecies beta\n"
            )
            hits = parse_speciesfinder_results(result)
        self.assertEqual([hit.template_accession for hit in hits], ["R2"])

    def test_placeholder_taxa_cannot_be_species_candidates(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = Path(temporary) / "results.txt"
            result.write_text(
                "# Assembly\tScore\tQuery_Coverage\tTemplate_Coverage\t"
                "Accession Number\tDescription\tTAXID\tTaxonomy\tTAXID\tSpecies\n"
                "A\t1000\t100\t99\tR1\tEnterobacter isolate\t"
                "1\tBacteria\t1\tEnterobacter sp\n"
                "B\t900\t100\t98\tR2\tfamily isolate\t"
                "2\tBacteria\t2\tEnterobacteriaceae bacterium\n"
                "C\t800\t100\t97\tR3\tKosakonia isolate\t"
                "3\tBacteria\t3\tKosakonia oryzae\n"
            )
            hits = parse_speciesfinder_results(result)
        self.assertEqual([hit.template_accession for hit in hits], ["R3"])
        self.assertTrue(is_reportable_species_name("Candidatus Liberibacter asiaticus"))
        self.assertFalse(is_reportable_species_name("Uncultured bacterium"))

    def test_preliminary_signals_use_query_counts_when_available(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = Path(temporary) / "results.txt"
            result.write_text(
                "# Assembly\tNum\tScore\tQuery_Coverage\tTemplate_Coverage\t"
                "Accession Number\tDescription\tTAXID\tTaxonomy\t"
                "TAXID Species\tSpecies\n"
                "A\t900\t100\t99\t99\tR1\ta\t1\tBacteria\t1\tSpecies alpha\n"
                "B\t100\t50\t50\t50\tR2\tb\t2\tBacteria\t2\tSpecies beta\n"
            )
            hits = parse_speciesfinder_results(result)
        signals = speciesfinder_hits_to_signals(hits)
        self.assertEqual(signals[0].read_count, 900)
        self.assertEqual(signals[0].read_percent, 90.0)
        self.assertEqual(signals[0].nomination_score, 100.0)
        self.assertAlmostEqual(signals[0].nomination_percent, 100 * 100 / 150)
        self.assertIsNone(signals[0].distinct_minimizers)

    def test_real_kma_scores_are_not_mislabeled_as_read_counts(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = Path(temporary) / "results.txt"
            result.write_text(
                "# Assembly\tScore\tQuery_Coverage\tTemplate_Coverage\t"
                "Accession Number\tDescription\tTAXID\tTaxonomy\tTAXID\tSpecies\n"
                "A\t900\t99\t99\tR1\ta\t1\tBacteria\t1\tSpecies alpha\n"
                "B\t100\t50\t50\tR2\tb\t2\tBacteria\t2\tSpecies beta\n"
            )
            signals = speciesfinder_hits_to_signals(
                parse_speciesfinder_results(result)
            )
        self.assertEqual(signals[0].read_count, 0)
        self.assertEqual(signals[0].read_percent, 0.0)
        self.assertEqual(signals[0].nomination_score, 900.0)
        self.assertEqual(signals[0].nomination_percent, 90.0)
        self.assertIsNone(signals[0].reference_breadth)

    def test_runner_rejects_silent_internal_kma_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_path = root / "reads.fastq"
            input_path.write_text("@r\nACGT\n+\nIIII\n")
            prefix = root / "bacteria"
            for suffix in (".comp.b", ".seq.b", ".length.b", ".name"):
                Path(f"{prefix}{suffix}").write_text("x")
            taxonomy = root / "bacteria.tax"
            taxonomy.write_text("x")

            def fake_run(_command, **_kwargs):
                output = root / "output"
                output.mkdir(exist_ok=True)
                (output / "speciesfinder.log").write_text(
                    "KMA failed with return code 1\n"
                )
                (output / "results.txt").write_text("header\n")
                (output / "results.res").write_text("header\n")
                (output / "output.json").write_text(json.dumps({}))
                return SimpleNamespace(returncode=0, stdout="", stderr="")

            with patch(
                "bactipipe.speciesfinder_evidence.subprocess.run",
                side_effect=fake_run,
            ):
                with self.assertRaisesRegex(RuntimeError, "internal KMA failure"):
                    run_speciesfinder(
                        inputs=[input_path],
                        output_directory=root / "output",
                        database_prefix=prefix,
                        taxonomy_file=taxonomy,
                    )

    def test_runner_does_not_accept_stale_results_from_a_failed_retry(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_path = root / "reads.fastq"
            input_path.write_text("@r\nACGT\n+\nIIII\n")
            prefix = root / "bacteria"
            for suffix in (".comp.b", ".seq.b", ".length.b", ".name"):
                Path(f"{prefix}{suffix}").write_text("x")
            taxonomy = root / "bacteria.tax"
            taxonomy.write_text("x")
            output = root / "output"
            output.mkdir()
            (output / "results.txt").write_text("stale\n")
            (output / "results.res").write_text("stale\n")
            (output / "output.json").write_text("{}\n")

            with patch(
                "bactipipe.speciesfinder_evidence.subprocess.run",
                return_value=SimpleNamespace(returncode=0, stdout="", stderr=""),
            ):
                with self.assertRaisesRegex(RuntimeError, "required result"):
                    run_speciesfinder(
                        inputs=[input_path],
                        output_directory=output,
                        database_prefix=prefix,
                        taxonomy_file=taxonomy,
                    )

            self.assertFalse((output / "results.txt").exists())
            self.assertFalse((output / "results.res").exists())
            self.assertFalse((output / "output.json").exists())


if __name__ == "__main__":
    unittest.main()
