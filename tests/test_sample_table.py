import unittest
from bactipipe.sample_table import is_footer,valid_accession


class SampleTableTests(unittest.TestCase):
    def test_controlled_footer_markers(self):
        for name in ('*Analysis type: R = Relatedness', 'Native Barcoding Kit V14:',
                     'Flow Cell:', 'BSA:', 'MinKnow Software Version:'):
            self.assertTrue(is_footer(name))
        self.assertTrue(is_footer('', 'Lot # (Kit)'))

    def test_sample_names_and_missing_species_are_not_footers(self):
        for sample in ('Lambda','P2637773-1','BSA','Control','S01'):
            self.assertFalse(is_footer(sample,''))
            self.assertTrue(valid_accession(sample))
        self.assertFalse(valid_accession('Worksheet instructions: mix the reagents'))


if __name__=='__main__':unittest.main()
