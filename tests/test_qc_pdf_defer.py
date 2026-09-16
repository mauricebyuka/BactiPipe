import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from bactipipe.scripts import process_data


class DeferPDFTests(unittest.TestCase):
    def test_defer_writes_identical_qc_table_without_rendering_pdf(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / 'temp.tsv'
            source.write_text('S1\t25\tPass\tUnknown\tNo confident identification\tN/A\t100\t50\tPass\tReview\n')
            tables = []
            for defer in (False, True):
                target = root / f'{defer}.tsv'
                with patch.object(process_data, 'render_pdf_from_tsv') as renderer, \
                     patch.object(process_data, 'simple_print'), patch.object(process_data, 'logger'):
                    process_data.make_summary(str(target), str(source), ['Run name: test'],
                                              {'S1': (99, 0)}, defer_pdf=defer)
                    self.assertEqual(renderer.call_count, 0 if defer else 1)
                tables.append(target.read_text())
            self.assertEqual(*tables)


if __name__ == '__main__':
    unittest.main()
