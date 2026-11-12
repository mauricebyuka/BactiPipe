import ast
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Barrier
from types import SimpleNamespace
import unittest
from unittest.mock import patch, Mock

from bactipipe.identification_parallel import identify_batch, result_for, plan_identification


class IdentificationParallelTests(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / 'database_identity.json').write_text(json.dumps({
            'database_kind': 'speciesfinder_kma', 'components': {'db': {'size_bytes': 1024}}}))

    @patch('bactipipe.identification_parallel.cpu_budget', side_effect=lambda n: n)
    def test_cpu_memory_and_unknown_database_limits(self, _):
        for cpus, memory, expected in ((60, 96, 2), (4, 96, 1), (60, 40, 1)):
            plan = plan_identification(9, cpus, self.root, memory_gb=memory)
            self.assertEqual(plan.jobs, expected)
            self.assertLessEqual(plan.jobs * plan.threads_per_sample, cpus)
        self.assertEqual(plan_identification(1, 60, self.root, memory_gb=96).jobs, 1)
        (self.root / 'database_identity.json').unlink()
        self.assertEqual(plan_identification(9, 60, self.root, memory_gb=96).jobs, 1)

    @patch('bactipipe.identification_parallel.cpu_budget', side_effect=lambda n: n)
    @patch('bactipipe.identification_parallel.emit')
    def test_jobs_overlap_fail_independently_and_keep_sample_keys(self, emit, _):
        gate = Barrier(2)
        def identify(sample, expected, genome, reads, threads):
            self.assertEqual(threads, 30)
            gate.wait(timeout=3)
            if sample == 'B':
                raise RuntimeError('Failed nomination')
            return 'A result'
        results = identify_batch({'A': ('Unknown', 'a.fa', []), 'B': ('Unknown', 'b.fa', [])},
            identify, cpus=60, database_root=self.root, memory_gb=96, report=lambda _: None,
            audit_path=self.root / 'resources.json')
        self.assertEqual(result_for(results, 'A'), 'A result')
        with self.assertRaisesRegex(RuntimeError, 'Failed nomination'):
            result_for(results, 'B')
        emit.assert_any_call('identification', 'completed', sample='A')
        emit.assert_any_call('identification', 'failed', sample='B')

    def test_both_cli_control_branches_emit_terminal_progress(self):
        # Execute the actual Lambda branches without invoking full QC tools.
        root = Path(__file__).parents[1] / 'bactipipe' / 'scripts'
        for platform in ('nanopore', 'illumina'):
            tree = ast.parse((root / f'run_{platform}_qc.py').read_text())
            branch = next(n for n in ast.walk(tree) if isinstance(n, ast.If)
                and isinstance(n.test, ast.Call) and isinstance(n.test.func, ast.Name)
                and n.test.func.id == 'is_lambda_control')
            for status, terminal in [('Pass', 'completed'), ('Fail', 'completed'), ('Error', 'failed')]:
                with self.subTest(platform=platform, status=status):
                    emit = Mock()
                    env = dict(is_lambda_control=lambda _: True, organism='Lambda', sample='L', genome='L.fa',
                        qc_out=str(self.root), os=__import__('os'), args=SimpleNamespace(lambda_min_identity=99, lambda_min_coverage=90),
                        assess_control=lambda *a, **k: {'status': status, 'reason': 'test'},
                        writer=Mock(), control_summary_row=Mock(), avqc=20, qc_verdict='Pass',
                        cov_display='100X', mincov=30, min_cov=30, cov_verdict='Pass',
                        sample_failure_audit=Mock(), sample_failures={}, reported_samples=set(), emit=emit)
                    loop = ast.For(target=ast.Name(id='_', ctx=ast.Store()), iter=ast.List(elts=[ast.Constant(1)], ctx=ast.Load()), body=[branch], orelse=[])
                    exec(compile(ast.fix_missing_locations(ast.Module(body=[loop], type_ignores=[])), '<control>', 'exec'), env)
                    emit.assert_called_once_with('identification', terminal, sample='L')
                    self.assertEqual(env['reported_samples'], {'L'})

    def test_illumina_low_quality_branch_emits_skipped(self):
        source = Path(__file__).parents[1] / 'bactipipe/scripts/run_illumina_qc.py'
        branch = next(n for n in ast.walk(ast.parse(source.read_text()))
            if isinstance(n, ast.If) and isinstance(n.test, ast.Compare)
            and isinstance(n.test.left, ast.Name) and n.test.left.id == 'qc_verdict')
        emit = Mock()
        env = dict(sample='S', avqc=10, qc_verdict='Fail', organism='Unknown', min_cov=50,
                   reported_samples=set(), writer=Mock(), emit=emit)
        exec(compile(ast.fix_missing_locations(ast.Module(body=branch.orelse, type_ignores=[])), '<low-quality>', 'exec'), env)
        emit.assert_called_once_with('identification', 'skipped', sample='S')
