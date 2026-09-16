from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from bactipipe.reporting_groups import ReportingGroupCatalog


class ReportingGroupTests(unittest.TestCase):
    def test_packaged_catalog_contains_ecoli_shigella_complex(self):
        catalog = ReportingGroupCatalog.from_package()
        self.assertEqual(
            catalog.display_name_for_taxid("562"),
            "Escherichia coli/Shigella complex",
        )
        self.assertIn(
            "Shigella sonnei",
            catalog.member_names_by_display_name[
                "Escherichia coli/Shigella complex"
            ],
        )

    def test_species_cannot_belong_to_two_groups(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "groups.tsv"
            path.write_text(
                "group_id\tdisplay_name\tspecies_taxid\tspecies_name\n"
                "one\tGroup one\t1\tSpecies alpha\n"
                "two\tGroup two\t1\tSpecies alpha\n"
            )
            with self.assertRaisesRegex(ValueError, "more than one reporting group"):
                ReportingGroupCatalog.from_path(path)


if __name__ == "__main__":
    unittest.main()
