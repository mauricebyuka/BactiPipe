from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import tempfile
import unittest
import zipfile
from pathlib import Path

from bactipipe.reference_sequences import (
    download_reference_genomes,
    write_reviewed_manifest,
)


class ReferenceSequenceTests(unittest.TestCase):
    def test_download_is_checkpointed_and_manifest_preserves_multicontig_fasta(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inventory = root / "inventory.tsv"
            write_tsv(
                inventory,
                [
                    inventory_row("GCF_000001.1", "1", "Species alpha", "canonical"),
                    inventory_row("GCF_000002.1", "1", "Species alpha", "diversity"),
                ],
            )
            calls = []

            def runner(command, *, timeout):
                calls.append((command, timeout))
                archive = Path(command[command.index("--filename") + 1])
                accessions = Path(command[command.index("--inputfile") + 1]).read_text().splitlines()
                with zipfile.ZipFile(archive, "w") as output:
                    for accession in accessions:
                        output.writestr(
                            f"ncbi_dataset/data/{accession}/{accession}_ASM_test_genomic.fna",
                            f">contig1\nACGT\n>contig2\nNNAA\n",
                        )
                return subprocess.CompletedProcess(command, 0, "", "")

            output = root / "download"
            summary = download_reference_genomes(
                inventory,
                output,
                datasets_executable="/bin/true",
                batch_size=2,
                runner=runner,
            )
            self.assertEqual(summary["downloaded_count"], 2)
            self.assertEqual(summary["missing_count"], 0)
            self.assertEqual(len(calls), 1)
            audit = list(
                csv.DictReader(
                    (output / "sequence_download_audit.tsv").open(), delimiter="\t"
                )
            )
            self.assertEqual(audit[0]["sequence_count"], "2")
            self.assertEqual(audit[0]["total_bases"], "8")

            repeated = download_reference_genomes(
                inventory,
                output,
                datasets_executable="/bin/true",
                batch_size=2,
                runner=runner,
            )
            self.assertEqual(repeated["downloaded_count"], 2)
            self.assertEqual(len(calls), 1)

            manifest = output / "reviewed_reference_manifest.tsv"
            counts = write_reviewed_manifest(
                inventory, output / "sequence_download_audit.tsv", manifest
            )
            self.assertEqual(counts["manifest_references"], 2)
            rows = list(csv.DictReader(manifest.open(), delimiter="\t"))
            by_accession = {row["accession"]: row for row in rows}
            self.assertEqual(by_accession["GCF_000001.1"]["include_in_kraken"], "yes")
            self.assertEqual(by_accession["GCF_000002.1"]["include_in_kraken"], "no")
            self.assertEqual(by_accession["GCF_000001.1"]["reference_role"], "type")
            self.assertEqual(by_accession["GCF_000002.1"]["reference_role"], "diversity")
            self.assertTrue(
                (manifest.parent / by_accession["GCF_000001.1"]["fasta_path"]).is_file()
            )
            self.assertEqual(
                json.loads(manifest.with_suffix(".summary.json").read_text())[
                    "manifest_species"
                ],
                1,
            )

    def test_reviewed_manifest_applies_controlled_reporting_group(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inventory = root / "inventory.tsv"
            write_tsv(
                inventory,
                [inventory_row("GCF_000005845.2", "562", "Escherichia coli", "canonical")],
            )
            references = root / "references"
            references.mkdir()
            fasta = references / "GCF_000005845.2.fna.gz"
            fasta.write_bytes(b"placeholder")
            audit = root / "sequence_download_audit.tsv"
            write_tsv(
                audit,
                [
                    {
                        "accession": "GCF_000005845.2",
                        "status": "downloaded",
                        "fasta_path": "references/GCF_000005845.2.fna.gz",
                        "sequence_count": "1",
                        "total_bases": "8",
                        "sha256": hashlib.sha256(b"placeholder").hexdigest(),
                    }
                ],
            )
            manifest = root / "reviewed.tsv"
            summary = write_reviewed_manifest(inventory, audit, manifest)
            row = next(csv.DictReader(manifest.open(), delimiter="\t"))
            self.assertEqual(
                row["reporting_group"],
                "Escherichia coli/Shigella complex",
            )
            self.assertEqual(summary["reporting_group_members"], 1)

    def test_reviewed_manifest_preserves_unique_controlled_subspecies_name(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            row = inventory_row(
                "GCA_900637675.1", "1336", "Streptococcus equi", "canonical"
            )
            row["source_taxid"] = "148942"
            row["controlled_catalog_name"] = "Streptococcus equi subsp. equi"
            row["current_submitted_species"] = "Streptococcus equi"
            inventory = root / "inventory.tsv"
            write_tsv(inventory, [row])
            references = root / "references"
            references.mkdir()
            fasta = references / "GCA_900637675.1.fna.gz"
            fasta.write_bytes(b"placeholder")
            audit = root / "sequence_download_audit.tsv"
            write_tsv(
                audit,
                [
                    {
                        "accession": "GCA_900637675.1",
                        "status": "downloaded",
                        "fasta_path": "references/GCA_900637675.1.fna.gz",
                        "sequence_count": "1",
                        "total_bases": "8",
                        "sha256": hashlib.sha256(b"placeholder").hexdigest(),
                    }
                ],
            )
            manifest = root / "reviewed.tsv"
            write_reviewed_manifest(inventory, audit, manifest)
            manifest_row = next(csv.DictReader(manifest.open(), delimiter="\t"))
            self.assertEqual(
                manifest_row["species_name"],
                "Streptococcus equi subsp. equi",
            )


def inventory_row(accession: str, taxid: str, species: str, role: str) -> dict[str, str]:
    return {
        "accession": accession,
        "species_taxid": taxid,
        "source_species_taxid": taxid,
        "source_taxid": taxid,
        "ncbi_species_name": species,
        "gtdb_species_name": species,
        "source": "RefSeq",
        "assembly_level": "Complete Genome",
        "type_material": "True",
        "current_total_length": "8",
        "reference_role": role,
        "controlled_catalog_name": species,
        "controlled_match": "ncbi_species_taxid",
        "current_submitted_species": species,
    }


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
