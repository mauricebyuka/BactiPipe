"""Exercise log statements without importing the argument-parsing QC script."""
import ast
from pathlib import Path
import unittest
from unittest.mock import Mock


class IlluminaLogMessagesTests(unittest.TestCase):
    def messages(self, variable, **values):
        source = Path(__file__).resolve().parents[1] / 'bactipipe/scripts/run_illumina_qc.py'
        text = source.read_text()
        self.assertNotIn('Yield depth will be assessed after assembly.', text)
        results = []
        for node in ast.walk(ast.parse(text)):
            body = getattr(node, 'body', None)
            if not isinstance(body, list):
                continue
            for index, statement in enumerate(body):
                if not isinstance(statement, ast.Assign) or not any(
                    isinstance(target, ast.Name) and target.id == variable
                    for target in statement.targets
                ):
                    continue
                console, log = Mock(), Mock()
                scope = dict(values, simple_print=console, time_print=console,
                             logger=log, log='test.log')
                block = ast.Module(body=body[index:index + 3], type_ignores=[])
                exec(compile(block, str(source), 'exec'), scope)
                message = scope[variable]
                self.assertEqual(console.call_args.args[0], message)
                self.assertEqual(log.call_args.args[1], message)
                results.append(message)
        self.assertTrue(results)
        return results

    def test_quality_labels_sample_and_preserves_verdict(self):
        for verdict in ('Pass', 'Fail'):
            self.assertEqual(self.messages('quality_message', sample='S01', avqc=37.47,
                                           qc_verdict=verdict),
                             [f'Sample S01: mean quality score: {verdict} --- 37.47'])

    def test_measured_depth_and_final_required_depth(self):
        for coverage, verdict in (('150.25X', 'Pass'), ('N/A', 'N/A')):
            self.assertCountEqual(self.messages('depth_message', sample='S02',
                cov_display=coverage, cov_verdict=verdict, min_cov=40), [
                f'Sample S02: assembly-based yield depth: {coverage}',
                f'Sample S02: yield depth QC: {verdict} --- {coverage}; required 40X',
            ])


if __name__ == '__main__':
    unittest.main()
