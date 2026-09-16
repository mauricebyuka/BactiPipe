from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from bactipipe.ncbi_reference_metadata import (
    _resolve_batch_reports,
    fetch_assembly_reports,
    read_assembly_accessions,
)


class NcbiReferenceMetadataTests(unittest.TestCase):
    def test_invalid_accession_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "accessions.txt"
            source.write_text("NZ_NOT_AN_ASSEMBLY\n")
            with self.assertRaisesRegex(ValueError, "Invalid assembly accession"):
                read_assembly_accessions(source)

    def test_fetch_is_batched_audited_and_resumable(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            accessions = root / "accessions.txt"
            accessions.write_text("GCF_0001.1\nGCA_0002.1\nGCF_0001.1\n")
            counter = root / "calls.txt"
            executable = root / "datasets"
            executable.write_text(
                "#!/usr/bin/env python3\n"
                "import json, pathlib, sys\n"
                "if sys.argv[1] == 'version':\n"
                "    print('datasets version: test')\n"
                "    raise SystemExit(0)\n"
                f"counter = pathlib.Path({str(counter)!r})\n"
                "counter.write_text(counter.read_text() + 'call\\n' if counter.exists() else 'call\\n')\n"
                "source = pathlib.Path(sys.argv[sys.argv.index('--inputfile') + 1])\n"
                "for accession in source.read_text().splitlines():\n"
                "    if accession.startswith('GCF_'):\n"
                "        current = accession.rsplit('.', 1)[0] + '.2'\n"
                "        print(json.dumps({'accession': current, 'assembly_info': {'assembly_status': 'current'}}))\n"
            )
            executable.chmod(0o755)
            output = root / "reports"

            summary = fetch_assembly_reports(
                accessions,
                output,
                datasets_executable=str(executable),
                batch_size=1,
            )

            self.assertEqual(summary["requested_count"], 2)
            self.assertEqual(summary["retrieved_count"], 1)
            self.assertEqual(summary["missing_count"], 1)
            self.assertEqual(summary["current_version_substitution_count"], 1)
            self.assertEqual(summary["status"], "success_with_missing_accessions")
            self.assertEqual((output / "missing_accessions.txt").read_text(), "GCA_0002.1\n")
            self.assertEqual(len(counter.read_text().splitlines()), 2)
            self.assertTrue((output / "FETCH_SUCCESS").exists())

            repeated = fetch_assembly_reports(
                accessions,
                output,
                datasets_executable=str(executable),
                batch_size=1,
            )
            self.assertEqual(repeated["retrieved_count"], 1)
            self.assertEqual(len(counter.read_text().splitlines()), 2)
            combined = [json.loads(line) for line in (output / "assembly_data_report.jsonl").read_text().splitlines()]
            self.assertEqual([row["accession"] for row in combined], ["GCF_0001.2"])
            resolution = (output / "accession_resolution.tsv").read_text()
            self.assertIn("GCF_0001.1\tGCF_0001.2\tcurrent_version_substituted", resolution)

    def test_resolution_prefers_current_when_ncbi_returns_version_history(self):
        reports = {
            "GCA_0001.1": {
                "accession": "GCA_0001.1",
                "assembly_info": {"assembly_status": "previous"},
            },
            "GCA_0001.2": {
                "accession": "GCA_0001.2",
                "assembly_info": {"assembly_status": "current"},
            },
        }
        resolution = _resolve_batch_reports(
            ["GCA_0001.1"], reports, batch_root=Path("batch")
        )
        self.assertEqual(resolution, {"GCA_0001.1": "GCA_0001.2"})


if __name__ == "__main__":
    unittest.main()
