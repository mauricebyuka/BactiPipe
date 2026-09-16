import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from bactipipe.iterative_phylogeny import place_iteratively, atomic_json, FILES, valid_receipt


class IterativePlacementTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.assembly = self.root / "assembly.fasta"
        self.assembly.write_text(">contig\nAAAA\n")
        self.calls = []
        self.sequence = "ACD-EF"
        self.count = 119
        self.identity = {"release": "r232", "parameters": {}}
        self.args = dict(assembly=self.assembly, cache_root=self.root / "cache",
                         reference_root=self.root, metadata=self.root / "metadata.tsv")
        def run(command, log, env, heartbeat, *args):
            self.calls.append(command)
            heartbeat()
            output = Path(command[command.index("--out_dir") + 1])
            output.mkdir(parents=True, exist_ok=True)
            if "identify" in command:
                protein = output / "identify/intermediate_results/marker_genes/query/query_protein.faa"
                protein.parent.mkdir(parents=True, exist_ok=True)
                protein.write_text(">protein\nACDEF\n")
                (output / "placement.bac120.markers_summary.tsv").write_text(
                    "name\tnumber_unique_genes\tnumber_multiple_genes\nquery\t" + str(self.count) + "\t0\n")
            if "align" in command:
                (output / "placement.bac120.user_msa.fasta").write_text(">query\n" + self.sequence + "\n")
        def render(**kwargs):
            from bactipipe.phylogenetic_placement import PlacementResult
            output = kwargs["output_directory"]
            for name in FILES:
                (output / name).write_text("tree data")
            return PlacementResult(status="available", reason="supporting", marker_evidence={"classification": "DO NOT REUSE ANI"})
        for target, value in (("execute", run), ("render_bac120_placement", render)):
            patcher = patch("bactipipe.iterative_phylogeny." + target, side_effect=value)
            patcher.start()
            self.addCleanup(patcher.stop)
        patcher = patch("bactipipe.iterative_phylogeny.reference_identity", side_effect=lambda *args: json.loads(json.dumps(self.identity)))
        patcher.start(); self.addCleanup(patcher.stop)
        patcher = patch("bactipipe.iterative_phylogeny.tool_identity", return_value={"tools": "digest"})
        patcher.start(); self.addCleanup(patcher.stop)
        patcher = patch("bactipipe.iterative_phylogeny.subprocess.check_output", return_value="gtdbtk 2.7.2")
        patcher.start(); self.addCleanup(patcher.stop)

    def test_identical_assembly_skips_all_stages(self):
        first = place_iteratively(**self.args)
        self.assertEqual(first["reuse"], "none")
        self.assertEqual(len(self.calls), 3)
        self.assertIn("--genes", self.calls[-1])
        second = place_iteratively(**self.args)
        self.assertEqual(second["reuse"], "assembly")
        self.assertEqual(len(self.calls), 3)

    def test_identical_markers_skip_placement_but_refresh_evidence(self):
        place_iteratively(**self.args)
        self.assembly.write_text(">contig\nAAAAAA\n")
        self.count = 118
        result = place_iteratively(**self.args)
        self.assertEqual(result["reuse"], "markers")
        self.assertEqual(len(self.calls), 5)
        payload = json.loads((Path(result["directory"]) / "phylogenetic_placement.json").read_text())
        self.assertEqual(payload["marker_evidence"]["single_copy_markers"], 118)
        self.assertNotIn("classification", payload["marker_evidence"])
        self.assertEqual(list((self.root / "cache").glob("work-*")), [])

    def test_changed_markers_recompute_placement(self):
        first = place_iteratively(**self.args)
        self.assembly.write_text(">contig\nAAAATT\n")
        self.sequence = "ACD-EG"
        second = place_iteratively(**self.args)
        self.assertNotEqual(first["marker_key"], second["marker_key"])
        self.assertEqual(len(self.calls), 6)

    def test_changed_reference_invalidates_both_caches(self):
        place_iteratively(**self.args)
        self.identity["release"] = "r233"
        second = place_iteratively(**self.args)
        self.assertEqual(second["reuse"], "none")
        self.assertEqual(len(self.calls), 6)

    def test_changed_settings_invalidate_cache(self):
        place_iteratively(**self.args)
        place_iteratively(**self.args, threads=4)
        self.assertEqual(len(self.calls), 6)

    def test_changed_exact_ani_panel_invalidates_compact_tree_cache(self):
        first=place_iteratively(**self.args,preferred_accessions=['GCF_1.1'])
        second=place_iteratively(**self.args,preferred_accessions=['GCF_1.1'])
        self.assertEqual(second['reuse'],'assembly')
        third=place_iteratively(**self.args,preferred_accessions=['GCF_1.2'])
        self.assertNotEqual(first['marker_key'],third['marker_key'])
        self.assertEqual(len(self.calls),6)

    def test_corrupt_svg_not_reused(self):
        first = place_iteratively(**self.args)
        (Path(first["directory"]) / FILES[0]).write_text("broken")
        self.assertIsNone(valid_receipt(Path(first["directory"])))
        second = place_iteratively(**self.args)
        self.assertEqual(second["reuse"], "markers")
        self.assertIsNotNone(valid_receipt(Path(second["directory"])))

    def test_failed_stage_cleans_heavy_work_but_does_not_cache_success(self):
        with patch("bactipipe.iterative_phylogeny.execute", side_effect=RuntimeError("tool failed")):
            with self.assertRaises(RuntimeError):
                place_iteratively(**self.args)
        self.assertFalse((self.root / "cache/assemblies").exists())
        self.assertEqual(list((self.root / "cache").glob("work-*")), [])

    def test_missing_markers_do_not_run_classify_or_invent_counts(self):
        self.sequence = ""
        result = place_iteratively(**self.args)
        self.assertEqual(result["status"], "unavailable")
        self.assertEqual(len(self.calls), 2)


if __name__ == "__main__":
    unittest.main()
