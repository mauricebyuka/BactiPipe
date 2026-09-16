from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from bactipipe.mixed_signal_mapping import (
    _species_breadth,
    _species_mapping_metrics,
    _write_combined_reference,
)
from bactipipe.species_identification import ReferenceCatalog, TaxonomicSignal


class MixedSignalMappingTests(unittest.TestCase):
    def test_combined_reference_namespaces_contigs_and_coverage_by_species(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            refs = root / "refs"
            refs.mkdir()
            alpha = refs / "GCF_1.fna"
            beta = refs / "GCF_2.fna"
            alpha.write_text(">chromosome alpha\nAAAA\n>plasmid\nAA\n")
            beta.write_text(">chromosome beta\nCCCCCC\n")
            manifest = root / "reference_manifest.tsv"
            header = (
                "accession\tspecies_taxid\tspecies_name\tfasta_path\tsource\t"
                "reference_role\tassembly_level\tsequence_origin\ttype_material\t"
                "canonical_for_tree\tcallable_for_identification\treporting_group\t"
                "ncbi_taxonomy_name\tgtdb_taxonomy\tsha256\n"
            )
            manifest.write_text(
                header
                + self._row("GCF_1", "1", "Species alpha", alpha)
                + self._row("GCF_2", "2", "Species beta", beta)
            )
            catalog = ReferenceCatalog.from_tsv(manifest)
            selected = [
                (TaxonomicSignal("1", "Species alpha", 100, 60.0), catalog.records[0]),
                (TaxonomicSignal("2", "Species beta", 50, 30.0), catalog.records[1]),
            ]
            combined = root / "combined.fna"
            contigs = _write_combined_reference(
                combined, catalog=catalog, selected=selected
            )
            self.assertIn(">taxid_1__GCF_1__1", combined.read_text())

            coverage = root / "coverage.tsv"
            names = list(contigs)
            coverage.write_text(
                "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\t"
                "meandepth\tmeanbaseq\tmeanmapq\n"
                f"{names[0]}\t1\t4\t10\t4\t100\t1\t30\t60\n"
                f"{names[1]}\t1\t2\t0\t0\t0\t0\t0\t0\n"
                f"{names[2]}\t1\t6\t5\t3\t50\t1\t30\t60\n"
            )
            breadth = _species_breadth(coverage, contigs)
            self.assertAlmostEqual(breadth["1"], 100 * 4 / 6)
            self.assertEqual(breadth["2"], 50.0)
            metrics = _species_mapping_metrics(coverage, contigs)
            self.assertEqual(metrics["1"][1], 10)
            self.assertEqual(metrics["2"][1], 5)

    @staticmethod
    def _row(accession: str, taxid: str, species: str, fasta: Path) -> str:
        digest = hashlib.sha256(fasta.read_bytes()).hexdigest()
        return (
            f"{accession}\t{taxid}\t{species}\t{fasta}\tRefSeq\ttype\tcomplete\t"
            f"isolate\tyes\tyes\tyes\t\t{species}\td__Bacteria;s__{species}\t{digest}\n"
        )


if __name__ == "__main__":
    unittest.main()
