from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
import json
from types import SimpleNamespace
from unittest.mock import patch

from bactipipe.scripts.sample_update import run_qc_sample_update


class SampleUpdateTests(unittest.TestCase):
    def test_update_replaces_selected_outputs_and_summary_row(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            run_name = "run-1"
            target = root / run_name
            self._write_result(target, run_name, {"S1": "old-one", "S2": "old-two"})

            def fake_subprocess(command, check):
                outdir = Path(command[command.index("--outdir") + 1]) / run_name
                self._write_result(outdir, run_name, {"S2": "new-two"})
                return SimpleNamespace(returncode=0)

            with patch("bactipipe.scripts.sample_update.subprocess.run", side_effect=fake_subprocess):
                run_qc_sample_update(
                    script_path=root / "runner.py",
                    argv=["--run_name", run_name, "--outdir", str(root), "--only-sample", "S2", "--update-existing"],
                    outdir_parent=root,
                    run_name=run_name,
                    sample_names=["S2"],
                )

            self.assertEqual((target / "qc_out" / "S1" / "marker.txt").read_text(), "old-one")
            self.assertEqual((target / "qc_out" / "S2" / "marker.txt").read_text(), "new-two")
            summary = (target / f"{run_name}_qc_summary.tsv").read_text()
            self.assertIn("S1\told-one", summary)
            self.assertIn("S2\tnew-two", summary)
            temporary_summary = (target / "qc_out" / "temp_qc_summary.tsv").read_text()
            self.assertIn("S1\told-one", temporary_summary)
            self.assertIn("S2\tnew-two", temporary_summary)
            history = json.loads(
                (target / "sample_update_history.jsonl").read_text().splitlines()[-1]
            )
            self.assertEqual(history["samples"], ["S2"])
            audit = target / history["provenance"]
            self.assertIn(
                "S1\told-one", (audit / f"{run_name}_qc_summary.tsv").read_text()
            )
            self.assertNotIn(
                "S1\told-one",
                (
                    audit / "selected_reanalysis" / f"{run_name}_qc_summary.tsv"
                ).read_text(),
            )

    def test_failed_reanalysis_leaves_existing_result_unchanged(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            run_name = "run-1"
            target = root / run_name
            self._write_result(target, run_name, {"S1": "published"})
            with (
                patch(
                    "bactipipe.scripts.sample_update.subprocess.run",
                    return_value=SimpleNamespace(returncode=1),
                ),
                self.assertRaisesRegex(RuntimeError, "status 1"),
            ):
                run_qc_sample_update(
                    script_path=root / "runner.py",
                    argv=["--run_name", run_name, "--only-sample", "S1", "--update-existing"],
                    outdir_parent=root,
                    run_name=run_name,
                    sample_names=["S1"],
                )
            self.assertEqual((target / "qc_out" / "S1" / "marker.txt").read_text(), "published")
            self.assertFalse((target / "sample_update_history.jsonl").exists())

    @staticmethod
    def _write_result(root: Path, run_name: str, samples: dict[str, str]) -> None:
        for sample, value in samples.items():
            sample_dir = root / "qc_out" / sample
            sample_dir.mkdir(parents=True, exist_ok=True)
            (sample_dir / "marker.txt").write_text(value)
            genome = root / "assemblies" / "genomes" / f"{sample}.fasta"
            genome.parent.mkdir(parents=True, exist_ok=True)
            genome.write_text(f">{sample}\nACGT\n")
        (root / f"{run_name}_qc_summary.tsv").write_text(
            "BactiPipe>>v0.1.0\nSample ID\tMean Quality\n"
            + "".join(f"{sample}\t{value}\n" for sample, value in samples.items())
        )
        (root / "qc_out" / "temp_qc_summary.tsv").write_text(
            "Sample\tMean Quality\n"
            + "".join(f"{sample}\t{value}\n" for sample, value in samples.items())
        )
        (root / "sample_failures.tsv").write_text("sample_id\tstage\treason\n")


if __name__ == "__main__":
    unittest.main()
