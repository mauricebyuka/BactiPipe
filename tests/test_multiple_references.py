import csv
import json
import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from bactipipe.scripts import type_genomes as T, utils as U
from bactipipe.relatedness_columns import reference_column_groups, report_sections
from bactipipe.downstream_manifest import read_downstream_manifest, write_command_sample_sheet


class MultipleReferenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.log = logging.getLogger('multiple-reference-test')
        self.refs = [self.root / f'R{i}.fasta' for i in range(1, 4)]
        for path in self.refs:
            path.write_text(f'>{path.stem}\nACGTACGT\n')

    def test_cli_flags_and_reference_validation(self):
        args = T._build_argparser().parse_args([
            '--reference', str(self.refs[0]), '--reference-id', 'R1', '--reference-specimen', 'Vaccine A',
            '--reference', str(self.refs[1]), '--reference-id', 'R2', '--reference-specimen', 'Vaccine B'])
        specs = T.resolve_reference_arguments(args)
        self.assertEqual([x[1:] for x in specs], [('R1', 'Vaccine A'), ('R2', 'Vaccine B')])
        args.reference_id.pop()
        with self.assertRaisesRegex(ValueError, 'once per'): T.resolve_reference_arguments(args)
        args.reference_id = None; args.reference_specimen = None
        args.reference = [self.refs[0], self.refs[0]]
        with self.assertRaisesRegex(ValueError, 'unique'): T.resolve_reference_arguments(args)
        args.reference = [self.root/'missing.fasta']
        with self.assertRaises(FileNotFoundError): T.resolve_reference_arguments(args)

    def test_manifest_multiple_references_excluded_from_command_sheet(self):
        path = self.root/'manifest.tsv'
        path.write_text('cohort_id\tcommand\torganism\tsample_id\trole\treference_name\n'
                        'P1\trelate\tsalmonella\tR1\treference\tR1\n'
                        'P1\trelate\tsalmonella\tR2\treference\tR2\n'
                        'P1\trelate\tsalmonella\tS1\tsample\t\n')
        c = read_downstream_manifest(path)[0]
        self.assertEqual([r.sample_id for r in c.references], ['R1','R2'])
        text = write_command_sample_sheet(c, self.root/'command.tsv').read_text()
        self.assertIn('S1', text); self.assertNotIn('R1', text); self.assertNotIn('R2', text)

    def test_summary_distinct_values_and_pdf_continuations(self):
        manifest = [T.SampleRecord(sid, 'reference' if sid.startswith('R') else 'Isolate A',
                                  'Vaccine <A> & B' if sid.startswith('R') else 'brain',
                                  self.root/f'{sid}.fasta', True, '') for sid in ['R1','R2','R3','S1']]
        comparisons = []
        for i, ref in enumerate(self.refs):
            ani, ska = self.root/f'ani{i}.tsv', self.root/f'ska{i}.tsv'
            ani.write_text('sample\tANI\tAF\ttool\n' + ''.join(f'{rec.sample}\t{99-i}.25\t80\tskani\n' for rec in manifest))
            ska.write_text('sample\tsnps\n' + ''.join(f'{rec.sample}\t{1000+i}\n' for rec in manifest))
            comparisons.append(dict(sample_id=ref.stem, ani_tsv=ani, ska_vs_ref=ska))
        path = U.write_final_summary(accession='P1', manifest=manifest, paths={'outdir':self.root},
                                     logger=self.log, reference_comparisons=comparisons)
        rows = list(csv.reader(path.open(), delimiter='\t'))
        self.assertEqual(rows[1][-6:], ['R1: SNP distance','R1: ANI (%)','R2: SNP distance','R2: ANI (%)','R3: SNP distance','R3: ANI (%)'])
        self.assertEqual(rows[5][-6:], ['1,000','99.25','1,001','98.25','1,002','97.25'])
        self.assertEqual(rows[2][8:10], ['-','-'])
        self.assertEqual(rows[3][10:12], ['-','-'])
        self.assertEqual(reference_column_groups(rows[1]), (8,['R1','R2','R3']))
        sections = list(report_sections([('Summary', rows[1:6])]))
        self.assertEqual(len(sections), 2)
        self.assertEqual(len(sections[0][1][0]), 12)
        self.assertTrue(sections[1][2])
        self.assertEqual(sections[1][1][-1][:8], rows[5][:8])
        U.render_pdf_type_genomes_full(final_tsv=path, out_pdf=self.root/'report.pdf', title='Relatedness',
            header_text='Run: Test\nAccession: P1', tool_versions={'BactiPipe':'0.1.0'}, run_name='Test', accession='P1')
        self.assertGreater((self.root/'report.pdf').stat().st_size, 2000)

    def test_pipeline_runs_typing_and_ska_once_with_multiple_ani_comparisons(self):
        sheet = self.root/'samples.tsv'; sheet.write_text('S1\tIsolate A\tbrain\n')
        (self.root/'S1.fasta').write_text('>S1\nACGT\n')
        (self.root/'mlst_db').mkdir()
        out = self.root/'output'
        def ska(**kwargs):
            path = kwargs['ska_root']/'P1_ska.distances.tsv'
            path.parent.mkdir(exist_ok=True)
            path.write_text('Sample1\tSample2\tDistance\nR1\tR2\t3\nR1\tS1\t4\nR2\tS1\t7\n')
            return path, ''
        def ani(**kwargs):
            ref = kwargs['manifest'][kwargs['reference_index']].sample
            kwargs['out_tsv'].write_text('sample\tANI\tAF\ttool\nS1\t' + ('99.9' if ref=='R1' else '98.8')+'\t80\tskani\n')
        argv = ['--name','Test','--run-name','Test','--accession','P1','--organism','listeria',
                '--assemblies-dir',str(self.root),'--sample-sheet',str(sheet),'--outdir',str(out),
                '--db-dir',str(self.root),'--no-pdf']
        for ref in self.refs[:2]:
            argv += ['--reference',str(ref),'--reference-id',ref.stem,'--reference-specimen','Vaccine '+ref.stem]
        with patch.object(T, '_check_tool', return_value=True), \
             patch.object(T, 'run_mlst', return_value='1') as mlst, \
             patch.object(T, 'write_reference_ani_table', side_effect=ani) as ani_run, \
             patch.object(T.ska2, 'run_ska_allpairs', side_effect=ska) as ska_run:
            self.assertEqual(T.main(argv), 0)
        self.assertEqual(mlst.call_count, 3)
        self.assertEqual(ska_run.call_count, 1)
        self.assertEqual(ani_run.call_count, 2)
        evidence = json.loads((out/'analysis_result.json').read_text())
        self.assertEqual(evidence['reference_samples'], ['R1','R2'])
        self.assertEqual([s['is_reference'] for s in evidence['samples']], [True,True,False])
        summary = (out/'P1_relate.tsv').read_text()
        self.assertIn('99.90',summary); self.assertIn('98.80',summary)
        self.assertIn('Vaccine R2',summary)
        with (out/'P1_relate.tsv').open() as handle:
            rows = list(csv.reader(handle, delimiter='\t'))
        self.assertEqual(rows[4][-4:], ['4','99.90','7','98.80'])


if __name__ == '__main__': unittest.main()
