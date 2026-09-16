import json
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace
import tempfile
import threading
import unittest

from bactipipe.nanopore_parallel import plan_flye_resources, run_assembly_jobs


def fake_assembly(job):
    if job.sample == 'bad':
        raise RuntimeError('simulated assembly failure')
    output = Path(job.assembly_dir)/'genomes'/f'{job.sample}.fasta'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(f'>process_group_{os.getpgrp()}\nACGT\n')


class NanoporeParallelTests(unittest.TestCase):
    def test_default_six_by_ten_not_division_by_total_samples(self):
        for count in (6, 24, 72):
            result = plan_flye_resources(count, 60, memory_gb=90)
            self.assertEqual((result.jobs, result.threads_per_sample), (6, 10))

    def test_small_batches_explicit_eight_and_cpu_memory_caps(self):
        for count, cpus, memory, requested, expected in [
            (2, 60, 90, None, (2, 10)), (24, 60, 90, 8, (6, 8)),
            (24, 8, 90, None, (1, 8)), (24, 60, 32, None, (4, 10)),
            (24, 60, 12, None, (1, 10)), (24, 60, 90, 60, (1, 60)),
        ]:
            result = plan_flye_resources(count, cpus, requested_threads=requested, memory_gb=memory)
            self.assertEqual((result.jobs, result.threads_per_sample), expected)
            self.assertLessEqual(result.jobs * result.threads_per_sample, cpus)

    def test_invalid_resources(self):
        for kwargs in ({'samples': 0}, {'total_cpus': 0}, {'requested_threads': 0},
                       {'memory_gb': 3}, {'memory_gb': float('nan')}):
            with self.assertRaises(ValueError):
                plan_flye_resources(**dict(dict(samples=24, total_cpus=60, memory_gb=90), **kwargs))

    def test_slow_first_sample_does_not_hide_later_completion(self):
        with tempfile.TemporaryDirectory() as tmp:
            jobs = [SimpleNamespace(sample=name, assembly_dir=str(Path(tmp)/'assemblies')) for name in ('slow', 'fast')]
            released = threading.Event()
            messages = []
            def worker(job):
                if job.sample == 'slow':
                    if not released.wait(5): raise RuntimeError('Completion reporter blocked')
                fake_assembly(job)
            def report(message):
                messages.append(message)
                if 'sample fast:' in message: released.set()
            state = run_assembly_jobs(jobs, workers=2, worker=worker, report=report,
                state_path=Path(tmp)/'progress.json', executor_factory=ThreadPoolExecutor)
            self.assertIn('1/2 finished', messages[0])
            self.assertIn('sample fast: completed', messages[0])
            self.assertEqual(state['completed'], 2)
            self.assertEqual(json.loads((Path(tmp)/'progress.json').read_text()), state)

    def test_real_pool_isolates_failure_and_inherits_cancellation_process_group(self):
        with tempfile.TemporaryDirectory() as tmp:
            jobs = [SimpleNamespace(sample=name, assembly_dir=str(Path(tmp)/'assemblies')) for name in ('ok', 'bad')]
            state = run_assembly_jobs(jobs, workers=2, worker=fake_assembly, report=lambda message: None,
                state_path=Path(tmp)/'progress.json')
            self.assertEqual((state['finished'], state['completed'], state['failed']), (2, 1, 1))
            self.assertIn(str(os.getpgrp()), (Path(tmp)/'assemblies/genomes/ok.fasta').read_text())
            self.assertIn('simulated assembly failure', (Path(tmp)/'qc_out/bad/assembly_failure.txt').read_text())
