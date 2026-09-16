import csv
import json
import tempfile
import unittest
from pathlib import Path

from bactipipe.assembly_depth import assembly_depth
from bactipipe.scripts.sample_failures import publish_failures


class AssemblyDepthTests(unittest.TestCase):
    def test_depth_uses_all_contigs_and_retained_bases_without_taxonomy(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / 'sample.fasta'
            fasta.write_text('>chromosome\n' + 'A' * 400 + '\n>other\n' + 'N' * 100 + '\n')
            audit = root / 'qc_out/S1/depth_metrics.json'
            depth, verdict = assembly_depth(25000, fasta, 50, audit_path=audit)
            self.assertEqual((depth, verdict), (50.0, 'Pass'))
            record = json.loads(audit.read_text())
            self.assertEqual(record['assembly_length_bp'], 500)
            self.assertEqual(record['assembly_contigs'], 2)
            self.assertEqual(record['retained_bases'], 25000)
            self.assertEqual(assembly_depth(24999, fasta, 50)[1], 'Fail')
            self.assertEqual(assembly_depth(0, fasta, 50), (0.0, 'Fail'))
            for missing in [None, 'N/A', float('nan'), float('inf'), -1]:
                self.assertEqual(assembly_depth(missing, fasta, 50), (None, 'N/A'))

    def test_missing_empty_or_invalid_assembly_does_not_use_an_estimate(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / 'sample.fasta'
            self.assertEqual(assembly_depth(250000000, fasta, 50), (None, 'N/A'))
            for contents in ['', '>empty\n', 'ACTG\n']:
                fasta.write_text(contents)
                self.assertEqual(assembly_depth(250000000, fasta, 50), (None, 'N/A'))

    def test_identification_execution_failure_keeps_already_measured_depth(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / 'sample.fasta'
            fasta.write_text('>one\nAAAA\n')
            assembly_depth(200, fasta, 50, audit_path=root / 'qc_out/S1/depth_metrics.json')
            summary = root / 'summary.tsv'
            publish_failures(root, summary, {'S1': 'Unknown'},
                             {'S1': ('organism_identification', 'tool failed')}, required_depth=50)
            with summary.open() as handle:
                row = list(csv.reader(handle, delimiter='\t'))[0]
            self.assertEqual(row[4], 'Analysis failed')
            self.assertEqual(row[6:9], ['50.00X', '50.0', 'Pass'])

    def test_both_qc_paths_use_the_assembly_based_final_metric(self):
        import ast
        root = Path(__file__).resolve().parents[1] / 'bactipipe/scripts'
        for script in ['run_illumina_qc.py', 'run_nanopore_qc.py']:
            tree = ast.parse((root / script).read_text())
            calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                     and isinstance(n.func, ast.Name) and n.func.id == 'assembly_depth']
            self.assertTrue(calls)
            for call in calls:
                self.assertEqual(call.args[0].id, 'total_bases')
                self.assertEqual(call.args[1].id, 'genome')
