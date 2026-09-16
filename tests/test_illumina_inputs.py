import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

from bactipipe.illumina_inputs import stage_s3_reads


class IlluminaStagingTests(unittest.TestCase):
    def test_source_alias_and_checkpoint_are_independent_of_output_name(self):
        with tempfile.TemporaryDirectory() as directory:
            def download(**kwargs):
                self.assertEqual(kwargs['prefix'], 'ngslab/DATA/illuminaData/source')
                sample = Path(directory) / 'source' / 'S01'
                sample.mkdir(parents=True, exist_ok=True)
                (sample / 'S01_R1.fastq.gz').write_bytes(b'fixture')
            downloader = Mock(side_effect=download)
            args = dict(bucket='private', prefix='/ngslab/', run_name='X_source',
                        input_run_name='source', target_root=directory, download=downloader)
            # An unrelated directory must not count as a completed download.
            (Path(directory) / 'X_source').mkdir()
            self.assertEqual(stage_s3_reads(**args), str(Path(directory) / 'source'))
            stage_s3_reads(**args)
            self.assertEqual(downloader.call_count, 1)

    def test_empty_source_has_actionable_error_and_no_completion_marker(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, 's3://private/DATA/illuminaData/source/'):
                stage_s3_reads(bucket='private', prefix='', run_name='source',
                               input_run_name=None, target_root=directory, download=Mock())
            self.assertFalse(list(Path(directory).rglob('*.json')))

    def test_interrupted_download_is_not_reused(self):
        with tempfile.TemporaryDirectory() as directory:
            sample = Path(directory) / 'source' / 'S01'
            sample.mkdir(parents=True)
            (sample / 'S01_R1.fastq.gz').write_bytes(b'partial')
            downloader = Mock(side_effect=OSError('interrupted'))
            with self.assertRaisesRegex(OSError, 'interrupted'):
                stage_s3_reads(bucket='private', prefix='', run_name='X_source',
                               input_run_name='source', target_root=directory, download=downloader)
            downloader.assert_called_once()
            self.assertFalse(list(Path(directory).rglob('*.json')))

    def test_reject_path_as_alias(self):
        with self.assertRaises(ValueError):
            stage_s3_reads(bucket='private', prefix='', run_name='output',
                           input_run_name='../source', target_root='/unused', download=Mock())
