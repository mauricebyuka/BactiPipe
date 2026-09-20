"""Tool-routing regression tests; no databases or analysis jobs are modified."""
import contextlib
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from bactipipe.scripts import run_traits, traits_amr, traits_db, traits_vf
from bactipipe.scripts.traits_organisms import (
    amrfinder_organism, parse_amrfinder_organisms, virulence_databases,
)


class SerialPool:
    def __init__(self, **kwargs): pass
    def map(self, function, values): return [function(value) for value in values]
    def close(self): pass
    def join(self): pass


class DetectAuditTests(unittest.TestCase):
    def test_list_parsing_does_not_turn_errors_into_organisms(self):
        for text in ("", "Software directory: /tmp\nERROR: Database missing", "Failed to open database"):
            with self.subTest(text=text), self.assertRaises(RuntimeError):
                parse_amrfinder_organisms(text)
        self.assertEqual(parse_amrfinder_organisms("", "Available --organism options: Salmonella,\n Escherichia\namrfinder took 0 seconds to complete"), ["Escherichia", "Salmonella"])
        self.assertEqual(parse_amrfinder_organisms("Salmonella\nEscherichia\n"), ["Escherichia", "Salmonella"])

    def test_failed_query_exits_nonzero_not_fake_list(self):
        failed = subprocess.CompletedProcess([], 1, "", "missing database")
        with patch.object(run_traits, "_conda_run", return_value=failed), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as exc:
            run_traits.main(["--organism-list"])
        self.assertEqual(exc.exception.code, 2)

    def test_list_is_plain_readable_names_without_tools(self):
        text = run_traits._format_organism_list(["Salmonella", "Salmonella_enterica", "Enterococcus_faecium"])
        self.assertIn('Enterococcus', text.splitlines())
        self.assertIn('Enterococcus faecium', text.splitlines())
        self.assertEqual(text.splitlines().count('Salmonella'), 1)
        self.assertIn('Other', text.splitlines())
        for token in ('\t', '_', 'AMRFinder', 'ABRicate', 'VirulenceFinder', 'ncbi', 'vfdb'):
            self.assertNotIn(token, text)

    def test_simple_names_roundtrip_to_tool_tokens(self):
        from bactipipe.scripts.traits_organisms import display_organism_name
        supported = ['Salmonella_enterica', 'Enterococcus_faecium', 'Staphylococcus_aureus']
        for token in supported:
            self.assertEqual(amrfinder_organism(display_organism_name(token), supported), token)
        self.assertEqual(amrfinder_organism('Salmonella enterica', ['Salmonella']), 'Salmonella')
        self.assertIsNone(amrfinder_organism('Enterococcus', supported))

    def test_exact_amr_match_and_genus_specific_virulence(self):
        supported = ["Enterococcus_faecium", "Salmonella"]
        self.assertIsNone(amrfinder_organism("Enterococcus", supported))
        self.assertEqual(amrfinder_organism("enterococcus faecium", supported), "Enterococcus_faecium")
        self.assertEqual(virulence_databases("Enterococcus cecorum"), "virulence_ent")
        self.assertIsNone(virulence_databases("Staphylococcus epidermidis"))

    def test_main_records_routes_and_honors_modes(self):
        cases = [("Enterococcus", [], False, True, True),
                 ("Enterococcus faecium", [], True, False, True),
                 ("Aeromonas", [], False, True, False),
                 ("Salmonella", [], True, False, False),
                 ("Enterococcus", ["--virulence-only"], False, False, True),
                 ("Enterococcus", ["--amr-only"], False, True, False)]
        for organism, flags, amr, ncbi, vf in cases:
            with self.subTest(organism=organism, flags=flags), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                (root / "S1.fasta").write_text(">S1\nACGT\n")
                (root / "sheet.tsv").write_text("S1\tIsolate 1\tspecimen\n")
                with patch.object(run_traits, "_list_amrfinder_organisms", return_value=["Enterococcus_faecium", "Salmonella"]) as query, patch.object(run_traits.mp, "Pool", SerialPool), patch.object(traits_amr, "run_amr_for_sample", return_value=[]) as amr_call, patch.object(traits_vf, "run_vf_for_sample", return_value=[]) as vf_call, patch.object(traits_db, "ensure_and_collect_versions", return_value={}) as versions, contextlib.redirect_stdout(io.StringIO()):
                    run_traits.main(["--analyst", "test", "--genomes-dir", tmp, "--sample-sheet", str(root / "sheet.tsv"), "--run-name", "test", "--organism", organism, "--outdir", str(root / "out"), "--no-pdf"] + flags)
                params = json.loads((root / "out/analysis_result.json").read_text())["parameters"]
                self.assertEqual((params["amrfinder"], params["abricate_ncbi"], params["virulencefinder"]), (amr, ncbi, vf))
                self.assertEqual(versions.call_args.kwargs["want_amr"], amr)
                self.assertEqual(versions.call_args.kwargs["want_ncbi"], ncbi)
                if "--virulence-only" in flags:
                    query.assert_not_called()
                    amr_call.assert_not_called()
                if "--amr-only" in flags:
                    vf_call.assert_not_called()

    def test_ncbi_command_and_failed_output_not_reused(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = dict(sample="S1", fasta="S1.fasta", out_root=tmp, min_id=85, min_cov=55,
                        enable_amrfinder=False, enable_ncbi=True, enable_resfinder=False,
                        enable_card=False, threads=2, organism="Enterococcus")
            def fake_run(command, **kwargs):
                self.assertIn("ncbi", command)
                self.assertNotIn("-O", command)
                self.assertIn("85", command)
                kwargs["stdout"].write("#FILE\tGENE\t%IDENTITY\t%COVERAGE\nS1\tgeneA\t99\t100\n")
                return subprocess.CompletedProcess(command, 0, "", "")
            with patch.object(traits_amr.subprocess, "run", side_effect=fake_run):
                rows = traits_amr.run_amr_for_sample(**args)
            self.assertEqual(rows[0]["database"], "ncbi")
            with patch.object(traits_amr.subprocess, "run", return_value=subprocess.CompletedProcess([], 1, "", "failed")), self.assertRaises(subprocess.CalledProcessError):
                traits_amr.run_amr_for_sample(**args)

    def test_enterococcus_vf_prefix_and_thresholds(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "virulencefinder_db").mkdir()
            args = dict(sample="S1", fasta="input", out_root=tmp, min_id=85, min_cov=55,
                        enable_vfinder=True, enable_vfdb=False, enable_ecoli_vf=False,
                        threads=1, logger_fn=lambda *a, **k: None, db_root=tmp, organism="Enterococcus")
            with patch.object(traits_vf, "_run_virulencefinder", return_value="results.tsv") as tool, patch.object(traits_vf, "_parse_vfinder_results", return_value=[]):
                self.assertEqual(traits_vf.run_vf_for_sample(**args), [])
            self.assertEqual(tool.call_args.kwargs, dict(min_id=.85, min_cov=.55, databases="virulence_ent"))
            with patch.object(traits_vf, "_run_virulencefinder", return_value=None), self.assertRaises(RuntimeError):
                traits_vf.run_vf_for_sample(**args)

    def test_amr_only_fallback_records_abricate_and_ncbi_versions(self):
        with patch.object(traits_db, "_abricate_version", return_value="v1"), patch.object(traits_db, "_abricate_list_versions", return_value={"ncbi": "date"}), patch.object(traits_db, "_amrfinder_versions") as amr:
            versions = traits_db.ensure_and_collect_versions(want_amr=False, want_vf=False, want_ncbi=True)
        self.assertEqual(versions, {"tools": {"abricate": "v1"}, "databases": {"ncbi": "date"}})
        amr.assert_not_called()


if __name__ == "__main__": unittest.main()
