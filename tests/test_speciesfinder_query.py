import csv
import hashlib
import tempfile
import unittest
from pathlib import Path

from bactipipe.speciesfinder_query import FRAGMENT_LIMIT_BP, prepare_nomination_query


class NominationQueryTests(unittest.TestCase):
    def test_lossless_boundaries_and_determinism(self):
        for length in (1, 19999, 20000, 20001, 40000, 40001):
            with self.subTest(length=length), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                sequence = ("acGTNrySWkmbdhv" * (length // 14 + 1))[:length]
                source = root / "assembly.fasta"
                source.write_text(">duplicate header\n" + sequence[:100] + "\n" + sequence[100:]
                                  + "\n>duplicate header\nNN\n")
                original = source.read_bytes()
                evidence = prepare_nomination_query(source, root / "first")
                prepare_nomination_query(source, root / "second")
                query = (root / "first/query.fasta").read_text()
                sequences = [line for line in query.splitlines() if not line.startswith(">")]
                self.assertEqual("".join(sequences), sequence + "NN")
                self.assertTrue(all(0 < len(seq) <= FRAGMENT_LIMIT_BP for seq in sequences))
                self.assertEqual(evidence["fragment_count"], (length + 19999) // 20000 + 1)
                self.assertEqual(evidence["total_bases"], length + 2)
                self.assertEqual(evidence["source_sha256"], hashlib.sha256(original).hexdigest())
                self.assertEqual(evidence["query_sha256"], hashlib.sha256(query.encode()).hexdigest())
                self.assertEqual(evidence["concatenated_sequence_sha256"],
                                 hashlib.sha256((sequence + "NN").encode()).hexdigest())
                self.assertEqual(source.read_bytes(), original)
                for filename in ("query.fasta", "query.json", "fragments.tsv"):
                    self.assertEqual((root / "first" / filename).read_bytes(),
                                     (root / "second" / filename).read_bytes())
                with (root / "first/fragments.tsv").open() as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(len({row["fragment_id"] for row in rows}), len(rows))
                for row, seq in zip(rows, sequences):
                    start, stop = int(row["start_0based"]), int(row["end_exclusive"])
                    original_seq = sequence if row["contig_number"] == "1" else "NN"
                    self.assertEqual(seq, original_seq[start:stop])

    def test_rejects_malformed_without_publishing_success_audit(self):
        for text in ("", "ACGT\n", ">\nACGT\n", ">a\n>b\nACGT\n", ">a\nAC-GT\n"):
            with self.subTest(text=text), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                source = root / "assembly.fasta"
                source.write_text(text)
                with self.assertRaises(ValueError):
                    prepare_nomination_query(source, root / "output")
                self.assertFalse((root / "output/query.json").exists())
                self.assertEqual(source.read_text(), text)

    def test_rejects_invalid_limits_and_existing_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "assembly.fasta"
            source.write_text(">a\nACGT\n")
            for limit in (0, -1, 1.5, True):
                with self.assertRaises(ValueError):
                    prepare_nomination_query(source, root / "output", limit_bp=limit)
            with self.assertRaises(FileExistsError):
                prepare_nomination_query(source, root)
            self.assertEqual(source.read_text(), ">a\nACGT\n")


if __name__ == "__main__":
    unittest.main()
