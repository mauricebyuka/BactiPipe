from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from bactipipe.species_database import (
    _fasta_suffix,
    _write_kraken_library,
    build_identification_database,
    validate_database,
)
from bactipipe.species_identification import ReferenceCatalog


class SpeciesDatabaseTests(unittest.TestCase):
    def test_fasta_suffix_excludes_accession_version(self):
        self.assertEqual(_fasta_suffix(Path("GCF_001027105.1.fna")), ".fna")
        self.assertEqual(_fasta_suffix(Path("GCF_001027105.1.fasta.gz")), ".fasta.gz")

    def test_build_is_immutable_and_inactive(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest = self._manifest(root)
            taxonomy = root / "taxonomy"
            taxonomy.mkdir()
            (taxonomy / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n")
            (taxonomy / "names.dmp").write_text(
                "1\t|\troot\t|\t\t|\tscientific name\t|\n"
            )

            def fake_run(command, **kwargs):
                if command[1] == "sketch":
                    output = Path(command[command.index("-o") + 1])
                    if not output.is_absolute():
                        output = Path(kwargs["cwd"]) / output
                    output.mkdir()
                    return SimpleNamespace(returncode=0, stdout="built\n", stderr="")
                if "--build" in command:
                    kraken = Path(command[command.index("--db") + 1])
                    for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
                        (kraken / name).write_text("built\n")
                return SimpleNamespace(returncode=0, stdout="skani 0.3.0\n", stderr="")

            with patch("bactipipe.species_database.subprocess.run", side_effect=fake_run):
                output = build_identification_database(
                    reviewed_manifest=manifest,
                    output_root=root / "database",
                    ncbi_taxonomy_dir=taxonomy,
                    version="test-v1",
                    threads=2,
                )

            identity = validate_database(output, allow_candidate=True)
            self.assertEqual(identity["status"], "candidate")
            self.assertFalse(identity["active"])
            self.assertTrue((output / "BUILD_SUCCESS").exists())
            self.assertTrue((output / "checksums.sha256").exists())
            self.assertEqual(
                (output / "skani_input_files.txt").read_text(),
                "references/GCF_1.fna\n",
            )
            self.assertEqual(identity["reporting_group_count"], 0)
            self.assertTrue((output / "kraken2" / "hash.k2d").exists())
            self.assertTrue((output / "ncbi_taxonomy" / "nodes.dmp").exists())
            self.assertFalse((output / ".kraken-input-batches").exists())
            with self.assertRaisesRegex(ValueError, "not an active production"):
                validate_database(output)

    def test_callable_mag_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest = self._manifest(root, origin="MAG")
            catalog = ReferenceCatalog.from_tsv(manifest)
            with self.assertRaisesRegex(ValueError, "not isolate-derived"):
                catalog.validate_files()

    def test_kraken_library_uses_explicit_manifest_taxids(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest = self._manifest(root)
            catalog = ReferenceCatalog.from_tsv(manifest)
            output = root / "kraken.fna"
            _write_kraken_library(
                output,
                database_root=root,
                records=list(catalog.records),
            )
            self.assertIn(">GCF_1_1|kraken:taxid|1", output.read_text())

    @staticmethod
    def _manifest(root: Path, *, origin: str = "isolate") -> Path:
        fasta = root / "GCF_1.fna"
        fasta.write_text(">GCF_1\nACGTACGT\n")
        checksum = hashlib.sha256(fasta.read_bytes()).hexdigest()
        manifest = root / "reviewed.tsv"
        manifest.write_text(
            "accession\tspecies_taxid\tspecies_name\tfasta_path\tsource\t"
            "reference_role\tassembly_level\tsequence_origin\ttype_material\t"
            "canonical_for_tree\tcallable_for_identification\treporting_group\t"
            "ncbi_taxonomy_name\tgtdb_taxonomy\tsha256\n"
            f"GCF_1\t1\tSpecies alpha\t{fasta}\tRefSeq\ttype\tcomplete\t{origin}\t"
            f"yes\tyes\tyes\t\tSpecies alpha\td__Bacteria;s__Species alpha\t{checksum}\n"
        )
        return manifest


if __name__ == "__main__":
    unittest.main()
