from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from bactipipe.kraken_evidence import parse_kraken2_report, write_taxonomic_signals
from bactipipe.species_identification import parse_taxonomic_signals


class KrakenEvidenceTests(unittest.TestCase):
    def test_parses_species_from_minimizer_report(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            report = root / "report.tsv"
            report.write_text(
                "90.0\t900\t0\t10000\t5000\tG\t10\t  Mannheimia\n"
                "80.0\t800\t700\t9000\t4500\tS\t75985\t    Mannheimia haemolytica\n"
                "6.0\t60\t50\t800\t300\tS\t85404\t    Mannheimia varigena\n"
            )
            signals = parse_kraken2_report(report)
            self.assertEqual([item.species_taxid for item in signals], ["75985", "85404"])
            self.assertEqual(signals[1].distinct_minimizers, 300)

            compact = root / "signals.tsv"
            write_taxonomic_signals(signals, compact)
            restored = parse_taxonomic_signals(compact)
            self.assertEqual(restored[0].species_name, "Mannheimia haemolytica")
            self.assertIsNone(restored[0].reference_breadth)

    def test_rejects_standard_report_without_minimizer_columns(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            report = Path(temporary_dir) / "report.tsv"
            report.write_text("80.0\t800\t700\tS\t75985\tMannheimia haemolytica\n")
            with self.assertRaisesRegex(ValueError, "expected 8"):
                parse_kraken2_report(report)

    def test_reviewed_subspecies_uses_direct_counts_without_parent_overlap(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            report = Path(temporary_dir) / "report.tsv"
            report.write_text(
                "20.0\t200\t80\t1000\t500\tS\t1336\t    Streptococcus equi\n"
                "12.0\t120\t120\t700\t400\tS1\t40041\t      "
                "Streptococcus equi subsp. zooepidemicus\n"
            )
            signals = parse_kraken2_report(
                report,
                identification_names_by_taxid={
                    "1336": "Streptococcus equi subsp. equi",
                    "40041": "Streptococcus equi subsp. zooepidemicus",
                },
            )
            by_taxid = {item.species_taxid: item for item in signals}
            self.assertEqual(by_taxid["1336"].read_count, 80)
            self.assertEqual(by_taxid["1336"].read_percent, 8.0)
            self.assertEqual(
                by_taxid["1336"].species_name,
                "Streptococcus equi subsp. equi",
            )
            self.assertEqual(by_taxid["40041"].read_count, 120)
            self.assertEqual(by_taxid["40041"].read_percent, 12.0)


if __name__ == "__main__":
    unittest.main()
