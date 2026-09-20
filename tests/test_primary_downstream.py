import csv
from io import StringIO
from pathlib import Path
import tempfile
import unittest

from openpyxl import Workbook
from bactipipe.primary_downstream import primary_downstream_tsv, organism_token, report_order
from bactipipe.downstream_manifest import _parse_rows, DownstreamManifestError, write_command_sample_sheet


class PrimaryDownstreamTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / 'arbitrary.xlsx'

    def sheet(self, samples, reference='RefA', extra=()):
        w = Workbook(); s = w.active
        for i in range(4): s.cell(5+i, 5, f'R{i+1}:')
        s.cell(5, 6, reference)
        for col, val in enumerate(('Accession Number', 'Host and Source', 'Bacteria Species',
                                   'Isolate ID', 'R/D*', 'NB Index'), 1):
            s.cell(10, col, val)
        for n, row in enumerate(samples, 11):
            for col, val in enumerate(row, 1): s.cell(n, col, val)
        for n, row in enumerate(extra, 13+len(samples)):
            for col, val in enumerate(row, 1): s.cell(n, col, val)
        w.save(self.path); w.close()

    def cohorts(self):
        value = primary_downstream_tsv(self.path, detect_organisms=('Salmonella', 'Escherichia'))
        return _parse_rows(csv.reader(StringIO(value), delimiter='\t')) if value else []

    def test_multiple_references_both_commands_and_metadata(self):
        self.sheet([('P2633261-6', 'tissue', 'Salmonella enterica', 'Z', 'R1,D'),
                    ('P2633261-3B', 'brain', 'Salmonella enterica Typhimurium', 'A', 'R1,D')],
                   'RefA; RefB\nRefC')
        cohorts = self.cohorts()
        self.assertEqual(len(cohorts), 2)
        relate = [c for c in cohorts if c.command == 'relate']
        self.assertEqual(len(relate), 1)
        self.assertEqual({r.reference_name for r in relate[0].references}, {'RefA','RefB','RefC'})
        for c in cohorts:
            self.assertTrue(c.cohort_id.startswith('P2633261-'))
            samples = [s for s in c.samples if s.role == 'sample']
            self.assertEqual([s.display_name for s in samples], ['A','Z'])
            self.assertEqual(samples[0].specimen, 'brain')
            self.assertEqual(samples[0].collection_date, '')
        self.assertEqual(relate[0].organism, 'salmonella')
        detect = next(c for c in cohorts if c.command=='detect')
        self.assertTrue(detect.options['detect_amr']); self.assertTrue(detect.options['detect_virulence'])
        self.assertFalse(detect.options['extended_virulence'])

    def test_empty_none_refs_and_worksheet_order_without_names(self):
        for reference in ('', 'None'):
            self.sheet([('P1-9','','Escherichia coli','','R1'), ('P1-1','','Escherichia coli','','R1')], reference)
            c = self.cohorts()[0]
            self.assertIsNone(c.reference)
            self.assertEqual([s.sample_id for s in c.samples], ['P1-9', 'P1-1'])
            self.assertTrue(c.options['run_ska'])

    def test_empty_requests_and_protocol_boundary(self):
        self.sheet([('Lambda','','Lambda','--','--'), ('P1-1','','Escherichia coli','','')],
                   extra=[('notes','','fake','','D')])
        self.assertEqual(self.cohorts(), [])

    def test_optional_column_old_workbook(self):
        w=Workbook(); w.active.append(['Accession Number','Bacteria Species']); w.active.append(['S1','Unknown'])
        w.save(self.path); w.close()
        self.assertEqual(self.cohorts(), [])

    def test_invalid_requests_and_no_unsafe_names(self):
        for request in ('R5', 'DR1', 'bad'):
            self.sheet([('P1-1','','Escherichia coli','',''+request)])
            with self.assertRaisesRegex(DownstreamManifestError, 'Row 11'): self.cohorts()
        self.sheet([('../bad','','Escherichia coli','','D')])
        with self.assertRaisesRegex(DownstreamManifestError, 'Accession Number'): self.cohorts()

    def test_other_mapping_does_not_borrow_species_specific_scheme(self):
        self.assertEqual(organism_token('Enterococcus cecorum', 'relate'), ('Other','Enterococcus'))
        self.assertEqual(organism_token('Enterococcus cecorum', 'detect', ('Enterococcus_faecium',)),
                         ('Other','Enterococcus'))
        self.assertEqual(organism_token('Salmonella enterica Typhimurium','relate'), ('salmonella',''))

    def test_stable_ids_and_named_first_unnamed_sheet_order(self):
        rows=[('P1-9','','Escherichia coli','','D'),('P1-7','','Escherichia coli','Z','D'),
              ('P1-1','','Escherichia coli','','D'),('P1-4','','Escherichia coli','A','D')]
        self.sheet(rows)
        c=self.cohorts()[0]
        self.assertEqual([s.sample_id for s in c.samples], ['P1-4','P1-7','P1-9','P1-1'])
        self.sheet(rows+[('OTHER','','Salmonella enterica','B','D')])
        self.assertEqual(self.cohorts()[0].cohort_id, c.cohort_id)

    def test_relate_reference_stays_first(self):
        self.sheet([('S1','','Salmonella enterica','A','R1')])
        c=self.cohorts()[0]
        self.assertEqual(c.samples[0].role, 'reference')
        command_sheet=write_command_sample_sheet(c, Path(self.tmp.name)/'command.tsv')
        self.assertNotIn('RefA',command_sheet.read_text())

    def test_single_no_reference_can_type_without_pairwise_ska(self):
        self.sheet([('S1','','Salmonella enterica','','R2')])
        c=self.cohorts()[0]
        self.assertIsNone(c.reference)
        self.assertFalse(c.options['run_ska'])

    def test_reference_set_reordering_does_not_split_or_requeue_cohort(self):
        rows = [('P1-1','','Salmonella enterica','','R1'), ('P1-2','','Salmonella enterica','','R1')]
        self.sheet(rows, 'RefB; RefA; RefA')
        first = self.cohorts()[0]
        self.sheet(rows, 'RefA,RefB')
        second = self.cohorts()[0]
        self.assertEqual(first.cohort_id, second.cohort_id)
        self.assertEqual(first.signature, second.signature)
        self.assertEqual(len(second.references), 2)


class ReportOrderingTests(unittest.TestCase):
    def test_readers_keep_missing_name_order_and_reference_description(self):
        import logging
        from bactipipe.scripts import type_genomes, run_traits, traits_report
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); sheet=root/'samples.tsv'
            sheet.write_text('S9\t\tbrain\nS7\tZ\ttissue\nS1\t\tfeces\nS4\tA\tblood\n')
            for sid in ('S9','S7','S1','S4'):
                (root/f'{sid}.fasta').write_text('>x\nACGT\n')
            rows=type_genomes._read_sample_sheet(sheet,logging.getLogger(__name__))
            self.assertEqual([r[0] for r in rows], ['S4','S7','S9','S1'])
            detected=run_traits._read_sample_sheet(str(sheet),str(root))
            self.assertEqual([r[0] for r in detected], ['S4','S7','S9','S1'])
            self.assertEqual(traits_report._samples_in_run(dict.fromkeys(['Z','A'])), ['Z','A'])
            manifest=type_genomes._build_manifest(rows,root,logging.getLogger(__name__))
            reference=root/'ref.fasta'; reference.write_text('>x\nACGT\n')
            manifest=type_genomes._prepend_reference_to_manifest(manifest,reference,'REF',
                logging.getLogger(__name__),specimen='Vaccine reference')
            self.assertEqual(manifest[0].specimen,'Vaccine reference')
            self.assertEqual(manifest[0].isolate,'reference')

    def test_web_order_is_not_resorted_by_display_fallback(self):
        import logging
        from bactipipe.scripts import type_genomes, run_traits
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); sheet=root/'samples.tsv'
            sheet.write_text('S9\tZ\tbrain\nS1\tS1\tfeces\n')
            for sid in ('S9','S1'): (root/f'{sid}.fasta').write_text('>x\nACGT\n')
            self.assertEqual(type_genomes._read_sample_sheet(sheet,logging.getLogger(__name__),sample_order='sheet')[0][0],'S9')
            self.assertEqual(run_traits._read_sample_sheet(str(sheet),str(root),sample_order='sheet')[0][0],'S9')


if __name__ == '__main__': unittest.main()
