import json
import os
import tempfile
import threading
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from bactipipe.illumina_parallel import (
    SampleJob, available_cpus, available_memory_gb, cpu_budget, plan_resources, prepare_batch,
    prepare_sample, _launch,
)
from bactipipe.illumina_filtering import assemble_selected


class IlluminaParallelTests(unittest.TestCase):
    def job(self, root, sample='S', assembler='Skesa'):
        return SampleJob(sample, ('R1.fq', 'R2.fq'), str(root / 'qc'),
                         str(root / 'assembly'), assembler, 'v1', 28, 15, 24)

    def test_affinity_quota_and_requested_cpu_bounds(self):
        with patch('os.sched_getaffinity', return_value=set(range(60))), \
             patch.object(Path, 'read_text', return_value='1500000 100000'):
            self.assertEqual(available_cpus(), 15)
            self.assertEqual(cpu_budget(60), 15)
            self.assertEqual(cpu_budget(4), 4)
            with self.assertRaises(ValueError): cpu_budget(0)

    def test_resources_for_both_assemblers(self):
        for assembler in ('Skesa', 'Spades'):
            with self.subTest(assembler=assembler):
                normal = plan_resources(9, 60, memory_gb=96, assembler=assembler)
                self.assertEqual((normal.jobs, normal.threads_per_sample, normal.memory_per_sample_gb), (3, 20, 32))
                four = plan_resources(9, 60, maximum_jobs=4, memory_gb=96, assembler=assembler)
                self.assertEqual((four.jobs, four.threads_per_sample), (4, 15))
                one = plan_resources(9, 60, maximum_jobs=1, maximum_threads=60, memory_gb=96, assembler=assembler)
                self.assertEqual((one.jobs, one.threads_per_sample), (1, 60))
                small = plan_resources(9, 4, memory_gb=12, assembler=assembler)
                self.assertEqual((small.jobs, small.threads_per_sample, small.memory_per_sample_gb), (1, 4, 12))
                ram = plan_resources(9, 60, memory_gb=48, assembler=assembler)
                self.assertEqual(ram.jobs, 2)
                self.assertLessEqual(ram.jobs * ram.memory_per_sample_gb, ram.memory_budget_gb)
        unicycler = plan_resources(9, 60, memory_gb=96, assembler='Unicycler')
        self.assertEqual((unicycler.jobs, unicycler.threads_per_sample), (1, 60))

    def test_memory_respects_worker_share_and_cgroup_headroom(self):
        gib = 1024 ** 3
        values = {'memory.max': str(32 * gib), 'memory.current': str(12 * gib)}
        with patch('bactipipe.illumina_parallel.available_cpus', return_value=60), \
             patch('os.cpu_count', return_value=120), \
             patch('psutil.virtual_memory', return_value=SimpleNamespace(available=200 * gib)), \
             patch.object(Path, 'read_text', autospec=True, side_effect=lambda path: values[path.name]):
            self.assertEqual(available_memory_gb(), 16)
            values['memory.max'] = 'max'
            self.assertEqual(available_memory_gb(), 80)

    def test_invalid_resources(self):
        for memory in (0, 3, float('nan'), float('inf')):
            with self.assertRaises(ValueError): plan_resources(2, 8, memory_gb=memory)
        for option in ({'maximum_jobs': 0}, {'maximum_threads': 0}):
            with self.assertRaises(ValueError): plan_resources(2, 8, memory_gb=96, **option)

    def test_bounded_concurrency_order_and_failed_sample(self):
        resources = plan_resources(6, 16, maximum_jobs=2, memory_gb=48)
        jobs = [self.job(Path('/unused'), str(i)) for i in range(6)]
        barrier = threading.Barrier(2, timeout=5)
        lock = threading.Lock()
        active = peak = 0
        def launch(job):
            nonlocal active, peak
            with lock:
                active += 1
                peak = max(peak, active)
            barrier.wait()
            with lock: active -= 1
            if job.sample == '2': raise RuntimeError('simulated sample failure')
            return dict(sample=job.sample, failure=None, genome=job.sample)
        messages = []
        results = prepare_batch(jobs, resources, launch=launch, report=messages.append)
        self.assertTrue(any('2 concurrent samples x' in message and 'allocated CPUs' in message for message in messages))
        self.assertTrue(any('6/6 finished' in message for message in messages))
        self.assertEqual(list(results), list(map(str, range(6))))
        self.assertEqual(peak, 2)
        self.assertIn('simulated sample failure', results['2']['failure']['reason'])
        self.assertEqual(sum(r['failure'] is None for r in results.values()), 5)
        with self.assertRaises(ValueError): prepare_batch([jobs[0], jobs[0]], resources, launch=launch)

    def test_worker_passes_per_sample_resources_for_both_assemblers(self):
        selected = {'selected': {'mean_quality': 31}}
        with tempfile.TemporaryDirectory() as tmp:
            for assembler in ('Skesa', 'Spades'):
                job = self.job(Path(tmp), assembler=assembler)
                with patch('bactipipe.illumina_filtering.select_reads', return_value=selected) as select, \
                     patch('bactipipe.illumina_filtering.assemble_selected', return_value='S.fasta') as assemble:
                    result = prepare_sample(job)
                self.assertIsNone(result['failure'])
                self.assertEqual(result['genome'], 'S.fasta')
                self.assertEqual(select.call_args.args[0], job.reads)
                self.assertEqual(assemble.call_args.args[5], 15)
                self.assertEqual(assemble.call_args.kwargs['memory_gb'], 24)
                state = json.loads((Path(job.qc_root) / 'S/illumina_worker_state.json').read_text())
                self.assertEqual(state['status'], 'completed')

    def test_low_quality_and_exceptions_are_sample_local(self):
        with tempfile.TemporaryDirectory() as tmp:
            job = self.job(Path(tmp))
            selected = {'selected': {'mean_quality': 27}}
            with patch('bactipipe.illumina_filtering.select_reads', return_value=selected), \
                 patch('bactipipe.illumina_filtering.assemble_selected') as assemble:
                result = prepare_sample(job)
                assemble.assert_not_called()
                self.assertIsNone(result['failure'])  # Measured QC failure, not execution failure.
                self.assertIsNone(result['genome'])
            selected['selected']['mean_quality'] = 31
            with patch('bactipipe.illumina_filtering.select_reads', return_value=selected), \
                 patch('bactipipe.illumina_filtering.assemble_selected', side_effect=RuntimeError('assembler failed')):
                result = prepare_sample(job)
                self.assertEqual(result['failure']['stage'], 'assembly')
                self.assertEqual(result['selected'], selected)
            with patch('bactipipe.illumina_filtering.select_reads', side_effect=ValueError('bad reads')):
                self.assertEqual(prepare_sample(job)['failure']['stage'], 'qc')

    def test_stale_result_rejected_and_process_group_inherited(self):
        with tempfile.TemporaryDirectory() as tmp:
            job = self.job(Path(tmp))
            root = Path(job.qc_root) / job.sample
            root.mkdir(parents=True)
            result_path = root / 'illumina_preparation.json'
            result_path.write_text('{"sample":"S", "failure":null}')
            with patch('subprocess.Popen', return_value=SimpleNamespace(returncode=0, stdout=[], wait=lambda: 0)) as run:
                with self.assertRaises(RuntimeError): _launch(job)
                self.assertNotIn('start_new_session', run.call_args.kwargs)
                self.assertNotIn('process_group', run.call_args.kwargs)
                self.assertEqual(run.call_args.kwargs['env']['OPENBLAS_NUM_THREADS'], '1')
                self.assertEqual(run.call_args.kwargs['env']['OMP_NUM_THREADS'], '1')
            self.assertFalse(result_path.exists())

    def test_assembler_cpu_memory_flags_single_and_paired(self):
        from bactipipe.scripts.process_data import assemble
        for assembler in ('Skesa', 'Spades'):
            for reads in (['R1.fq'], ['R1.fq', 'R2.fq']):
                with self.subTest(assembler=assembler, reads=reads), tempfile.TemporaryDirectory() as tmp:
                    def fake_run(cmd, **kwargs):
                        if '--version' in cmd:
                            return SimpleNamespace(returncode=0, stdout=b'v1', stderr=b'')
                        thread_flag = '--cores' if assembler == 'Skesa' else '--threads'
                        self.assertEqual(cmd[cmd.index(thread_flag) + 1], '15')
                        self.assertEqual(cmd[cmd.index('--memory') + 1], '24')
                        self.assertEqual(kwargs['env']['OMP_NUM_THREADS'], '15')
                        self.assertEqual(kwargs['env']['OPENBLAS_NUM_THREADS'], '1')
                        self.assertEqual(kwargs['env']['MKL_NUM_THREADS'], '1')
                        self.assertEqual(os.environ['OMP_NUM_THREADS'], '1')
                        self.assertNotIn('start_new_session', kwargs)
                        target = Path(tmp) / 'S' / ('assembly.fasta' if assembler == 'Skesa' else 'scaffolds.fasta')
                        target.parent.mkdir(exist_ok=True)
                        target.write_text('>contig\n' + 'A' * 500 + '\n')
                        return SimpleNamespace(returncode=0, stdout=b'assembler output\n')
                    with patch.dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'), \
                         patch('bactipipe.scripts.process_data.subprocess.run', side_effect=fake_run):
                        assemble('S', reads, tmp, assembler, 'Illumina', cpus=15, memory_gb=24, single=False)
                    self.assertTrue((Path(tmp) / 'genomes/S.fasta').is_file())
                    self.assertIn('assembler output', (Path(tmp) / 'S/S.log').read_text())

    def test_failed_assembly_log_retained_and_resource_cache_invalidation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            selection = dict(selection_id='same', selected={'total_bases': 100}, selected_paths=['R1.fq'])
            def fail(**kwargs):
                log = Path(kwargs['assembly_dir']) / 'S/S.log'
                log.parent.mkdir(parents=True)
                log.write_text('tool failed: out of memory')
                raise RuntimeError('tool failed')
            with patch('bactipipe.scripts.process_data.assemble', side_effect=fail):
                with self.assertRaises(RuntimeError):
                    assemble_selected('S', selection, root, 'Skesa', 'v1', 15, None, memory_gb=24)
            self.assertIn('out of memory', (root / 'S/S.log').read_text())
            self.assertFalse(list(root.glob('.illumina-*')))
            def success(**kwargs):
                genome = Path(kwargs['assembly_dir']) / 'genomes/S.fasta'
                genome.parent.mkdir(parents=True)
                genome.write_text('>S\nAAAA\n')
            with patch('bactipipe.scripts.process_data.assemble', side_effect=success) as assembler:
                for threads, memory in ((15, 24), (15, 24), (20, 24), (20, 32)):
                    assemble_selected('S', selection, root, 'Skesa', 'v1', threads, None, memory_gb=memory)
                self.assertEqual(assembler.call_count, 3)


if __name__ == '__main__':
    unittest.main()
