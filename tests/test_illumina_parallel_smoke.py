"""Opt-in real QC/assembler smoke; uses tiny generated reads, never run data.

BACTIPIPE_REAL_ASSEMBLERS=1 python -m unittest discover -s tests \
    -p test_illumina_parallel_smoke.py
"""
import json
import os
from pathlib import Path
import random
import tempfile
import unittest

from bactipipe.illumina_parallel import Resources, SampleJob, prepare_batch


@unittest.skipUnless(os.environ.get('BACTIPIPE_REAL_ASSEMBLERS') == '1', 'real-tool smoke is opt-in')
class ParallelSmoke(unittest.TestCase):
    def test_real_workers_both_assemblers_and_bad_sample(self):
        with tempfile.TemporaryDirectory(prefix='illumina-parallel-smoke-') as tmp:
            root = Path(tmp)
            rng = random.Random(1729)
            genome = ''.join(rng.choices('ACGT', k=20000))
            reads = root / 'input.fastq'
            with reads.open('w') as stream:
                for i in range(10000):
                    start = rng.randrange(len(genome) - 150)
                    sequence = genome[start:start + 150]
                    if i % 2:
                        sequence = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
                    stream.write(f'@read{i}\n{sequence}\n+\n' + 'I' * 150 + '\n')
            jobs = [SampleJob(assembler, (str(reads),), str(root / 'qc'), str(root / 'assemblies'),
                              assembler, 'smoke', 28, 2, 4) for assembler in ('Skesa', 'Spades')]
            jobs.append(SampleJob('bad', (str(root / 'missing.fastq'),), str(root / 'qc'),
                                  str(root / 'assemblies'), 'Skesa', 'smoke', 28, 1, 4))
            outcomes = prepare_batch(jobs, Resources(4, 2, 2, 4, 8), report=print)
            for job in jobs[:2]:
                result = outcomes[job.sample]
                self.assertIsNone(result['failure'], result)
                self.assertGreater(Path(result['genome']).stat().st_size, 19000)
                self.assertEqual(result['selected']['selected']['total_bases'], 1500000)
                state = Path(job.qc_root) / job.sample / 'illumina_worker_state.json'
                self.assertEqual(json.loads(state.read_text())['status'], 'completed')
            spades_log = (root / 'assemblies/Spades/Spades.log').read_text()
            self.assertIn('Maximum # of threads to use (adjusted due to OMP capabilities): 2', spades_log)
            self.assertNotIn('Maximum # of threads to use (adjusted due to OMP capabilities): 1', spades_log)
            print('Verified real SPAdes effective OpenMP threads: 2 (QC worker cap: 1)')
            self.assertEqual(outcomes['bad']['failure']['stage'], 'qc')
            self.assertEqual(list(outcomes), ['Skesa', 'Spades', 'bad'])


if __name__ == '__main__':
    unittest.main()
