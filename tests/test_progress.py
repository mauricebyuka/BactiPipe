import json
import unittest
from unittest.mock import patch
from bactipipe.progress import emit, PREFIX


class ProgressEmitterTests(unittest.TestCase):
    @patch.dict('os.environ', {'BACTIPIPE_PROGRESS_EVENTS': '1'})
    def test_event_is_atomic_and_bounded(self):
        with patch('os.write') as write:
            emit('assembly', 'completed', sample='x' * 900)
            payload = write.call_args.args[1]
            self.assertLess(len(payload), 4096)
            self.assertEqual(len(json.loads(payload.decode()[len(PREFIX):])['sample']), 200)
            write.assert_called_once()

    @patch.dict('os.environ', {'BACTIPIPE_PROGRESS_EVENTS': '1'})
    def test_disabled_and_broken_output_are_harmless(self):
        with patch('os.write', side_effect=OSError('pipe closed')) as write:
            emit('unknown')
            write.assert_not_called()
            emit('identification')

    @patch.dict('os.environ', {'BACTIPIPE_PROGRESS_EVENTS': '0'})
    def test_normal_cli_output_is_unchanged(self):
        with patch('os.write') as write:
            emit('identification')
            write.assert_not_called()
