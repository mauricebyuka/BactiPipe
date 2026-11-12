import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from bactipipe.lambda_control import is_lambda_control, control_summary_row, CONTROL_IDENTIFICATION, alignment_metrics, assess_control, percent_threshold
from bactipipe.qc_status import interpret
from bactipipe.scripts.process_data import make_summary, checkM_stats
from bactipipe.scripts.pdf_renderer import _build_table, LayoutConfig
import pandas as pd


class LambdaControlTests(unittest.TestCase):
    def test_report_survives_low_depth_bacteria_and_keeps_control_verdict(self):
        from bactipipe.batch_species import reconcile_batch_identifications
        from bactipipe.scripts.sample_failures import publish_failures
        for control_status in ('Pass', 'Fail'):
            with self.subTest(control_status=control_status), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                source, final = root / 'temp.tsv', root / 'summary.tsv'
                # One bacterium assembled below 50X; another could not assemble.
                rows = [
                    ['Sample', 'Q', 'qc', 'Expected', 'Hit', 'Match', 'Depth', 'Min', 'cov', 'tax'],
                    ['S1', '21', 'Pass', 'Salmonella enterica', 'Salmonella enterica', '99', '11.27X', '50', 'Fail', 'Pass'],
                    control_summary_row('Lambda', 21.9, 'Pass', 'Lambda', '300X', 50, 'Pass',
                        {'status': control_status, 'identity_percent': 99.996}),
                ]
                with source.open('w') as handle:
                    csv.writer(handle, delimiter='\t').writerows(rows)
                publish_failures(root, source,
                    {'Lambda': 'Lambda', 'S1': 'Salmonella enterica', 'S2': 'Salmonella enterica'},
                    {'S2': ('assembly', 'Genome assembly did not produce a usable FASTA file.')}, required_depth=50)
                # No species evidence directory must not prevent report publication.
                self.assertEqual(reconcile_batch_identifications(
                    species_output_root=root / 'species_identification', database_root=root / 'db',
                    expected_organisms={'Lambda': 'Lambda', 'S2': 'Salmonella enterica'},
                    platform='Nanopore'), {})
                make_summary(str(final), str(source), [], {'S1': ('99', '0')}, logfile=str(root / 'log'))
                lines = final.read_text().splitlines()
                start = next(i for i, line in enumerate(lines) if line.startswith('Sample ID\t'))
                results = list(csv.DictReader(lines[start:], delimiter='\t'))
                self.assertEqual(results[0]['Sample ID'], 'Lambda')
                self.assertEqual(results[0]['Overall Quality'], control_status)
                by_sample = {row['Sample ID']: row for row in results}
                self.assertEqual(by_sample['S1']['Coverage Depth'], '11.27X')
                self.assertEqual(by_sample['S1']['Required depth'], '50X')
                self.assertEqual(by_sample['S1']['Overall Quality'], 'Fail')
                self.assertEqual(by_sample['S2']['Overall Quality'], 'Fail')
                self.assertEqual(final.with_suffix('.pdf').read_bytes()[:5], b'%PDF-')

    def test_both_platforms_explicitly_assess_declared_controls(self):
        import ast
        for filename in ('run_illumina_qc.py', 'run_nanopore_qc.py'):
            source = Path(__file__).parents[1] / 'bactipipe' / 'scripts' / filename
            tree = ast.parse(source.read_text())
            branches = [node for node in ast.walk(tree) if isinstance(node, ast.If)
                        and isinstance(node.test, ast.Call) and isinstance(node.test.func, ast.Name)
                        and node.test.func.id == 'is_lambda_control']
            calls = {call.func.id for branch in branches for node in branch.body for call in ast.walk(node)
                     if isinstance(call, ast.Call) and isinstance(call.func, ast.Name)}
            self.assertIn('assess_control', calls, filename)
            self.assertIn('control_summary_row', calls, filename)

    def hsp(self, query, subject, start=1, end=None, score=100, name='contig'):
        end=end if end is not None else start+len(subject.replace('-',''))-1
        return f'{name}\tLambda\t1\t{len(query.replace("-",""))}\t{start}\t{end}\t0\t{score}\t{query}\t{subject}'

    def test_identity_and_coverage_are_separate(self):
        result=alignment_metrics([self.hsp('ACGT','ACGT')],10)
        self.assertEqual(result['identity_percent'],100)
        self.assertEqual(result['reference_coverage_percent'],40)
        empty=alignment_metrics([],10)
        self.assertIsNone(empty['identity_percent'])
        self.assertEqual(empty['reference_coverage_percent'],0)

    def test_overlaps_duplicates_and_segments_do_not_inflate_coverage(self):
        rows=[self.hsp('ACGT','ACGT'),self.hsp('ACGT','ACGT',name='duplicate'),
              self.hsp('GTAC','GTAC',start=3),self.hsp('AC','AC',start=9)]
        result=alignment_metrics(rows,10)
        self.assertEqual(result['covered_reference_bases'],8)
        self.assertEqual(result['reference_coverage_percent'],80)
        self.assertEqual(result['identity_percent'],100)
        self.assertEqual(alignment_metrics(rows[::-1],10),result)

    def test_reverse_strand_and_gaps(self):
        self.assertEqual(alignment_metrics([self.hsp('ACGT','ACGT',start=4,end=1)],4)['identity_percent'],100)
        deleted=alignment_metrics([self.hsp('AC-GT','ACTGT')],5)
        self.assertEqual(deleted['identity_percent'],80)
        self.assertEqual(deleted['reference_coverage_percent'],80)
        inserted=alignment_metrics([self.hsp('ACTGT','AC-GT')],4)
        self.assertEqual(inserted['identity_percent'],80)
        self.assertEqual(inserted['reference_coverage_percent'],100)
        ambiguous=alignment_metrics([self.hsp('ACNT','ACGT')],4)
        self.assertEqual(ambiguous['identity_percent'],75)
        self.assertEqual(ambiguous['reference_coverage_percent'],75)

    def test_bad_alignment_and_thresholds_fail_closed(self):
        with self.assertRaises(ValueError): alignment_metrics(['not a BLAST record'],10)
        for value in (0,101,float('nan'),float('inf')):
            with self.assertRaises(ValueError): percent_threshold(value)

    def test_explicit_control_evidence_allows_pass(self):
        self.assertEqual(interpret('Lambda','Lambda',['Pass','Pass','N/A','N/A'],control_qc='Pass')['overall'],'Pass')
        self.assertEqual(interpret('Lambda','Lambda not confirmed',['Pass']*4,control_qc='Fail')['overall'],'Fail')
        self.assertEqual(interpret('Lambda','Lambda',['Pass','N/A'],control_qc='Pass')['overall'],'Review')
        self.assertEqual(interpret('Lambda','Lambda',['Pass']*4,control_qc='Error')['overall'],'Fail')
        row=control_summary_row('L',30,'Pass','Lambda','100X',50,'Pass',{'status':'Pass','identity_percent':99.98})
        self.assertEqual(row[4:6],['Lambda','99.980000'])
        self.assertEqual(row[-1],'ControlPass')

    def test_tool_failure_retains_error_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); assembly=root/'assembly.fasta'; assembly.write_text('>seq\nACGT\n')
            with patch('bactipipe.lambda_control.subprocess.run',side_effect=FileNotFoundError('blastn missing')):
                result=assess_control(assembly,root/'control')
            self.assertEqual(result['status'],'Error')
            self.assertIsNone(result['identity_percent'])
            self.assertTrue((root/'control'/'result.json').is_file())
            self.assertTrue((root/'control'/'metrics.tsv').is_file())

    def test_control_requires_explicit_organism(self):
        self.assertTrue(is_lambda_control("  lAmBdA "))
        for organism in ("Unknown", "Escherichia coli", "Lambda-like", "", None):
            self.assertFalse(is_lambda_control(organism))

    def test_control_cannot_pass_without_confirmation(self):
        for hit in ("Lambda", "Escherichia coli", CONTROL_IDENTIFICATION):
            result=interpret("Lambda",hit,["Pass","Pass","Fail","Fail"],isolate_qc_failure=True)
            self.assertEqual(result['overall'],'Review')
            self.assertEqual(result['quality'],'Pass')
            self.assertEqual(result['identification'],'control_not_assessed')
        self.assertEqual(interpret('Lambda','Lambda',['Pass','Fail','Pass','Pass'])['overall'],'Fail')
        self.assertEqual(interpret('Lambda','Analysis failed',['N/A']*4)['overall'],'Fail')
        self.assertEqual(interpret('Lambda','Lambda',['Pass','N/A'])['quality'],'Not assessed')

    def test_summary_keeps_control_first_and_ignores_checkm(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); source=root/'temp.tsv'; final=root/'summary.tsv'
            rows=[['Sample','Q','qc','Expected','Hit','Match','Depth','Min','cov','tax'],
                  ['S1','30','Pass','Escherichia coli','Escherichia coli','99','100X','50','Pass','Pass'],
                  control_summary_row('Control01',30,'Pass','Lambda','100X',50,'Pass')]
            with source.open('w') as handle: csv.writer(handle,delimiter='\t').writerows(rows)
            with patch('bactipipe.scripts.process_data.render_pdf_from_tsv'), patch('bactipipe.scripts.process_data.simple_print'):
                make_summary(str(final),str(source),[],{'Control01':['0','100'],'S1':['99','0']},logfile=str(root/'log'))
            lines=final.read_text().splitlines(); start=next(i for i,l in enumerate(lines) if l.startswith('Sample ID\t'))
            result=list(csv.DictReader(lines[start:],delimiter='\t'))
            self.assertEqual([r['Sample ID'] for r in result],['Control01','S1'])
            control=result[0]
            self.assertEqual(control['Overall Quality'],'Review')
            self.assertEqual(control['Identified organism'],CONTROL_IDENTIFICATION)
            for metric in ('% Match','CheckM completeness','CheckM contamination'): self.assertEqual(control[metric],'N/A')
            self.assertEqual(control['Coverage Depth'],'100.00X')
            self.assertEqual(result[1]['Overall Quality'],'Pass')
            table=_build_table(pd.DataFrame(result),LayoutConfig())
            self.assertEqual(table._cellvalues[1][2].getPlainText(),'Lambda (Control)')
            self.assertEqual(table._cellvalues[1][7:9],['N/A','N/A'])

    def test_confirmed_and_partial_control_summary_exports(self):
        for status, expected_hit in [('Pass','Lambda'),('Fail','Lambda not confirmed')]:
            with self.subTest(status=status), tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp); source=root/'temp.tsv'; final=root/'summary.tsv'
                row=control_summary_row('Control01',30,'Pass','Lambda','100X',50,'Pass',
                                        {'status':status,'identity_percent':99.98})
                with source.open('w') as handle: csv.writer(handle,delimiter='\t').writerow(row)
                with patch('bactipipe.scripts.process_data.render_pdf_from_tsv'), patch('bactipipe.scripts.process_data.simple_print'):
                    make_summary(str(final),str(source),[],{},logfile=str(root/'log'))
                lines=final.read_text().splitlines(); start=next(i for i,l in enumerate(lines) if l.startswith('Sample ID\t'))
                result=list(csv.DictReader(lines[start:],delimiter='\t'))[0]
                self.assertEqual(result['Overall Quality'],status)
                self.assertEqual(result['Control QC'],status)
                self.assertEqual(result['% Match'],'99.98')
                self.assertEqual(result['Identified organism'],expected_hit)
                table=_build_table(pd.DataFrame([result]),LayoutConfig())
                self.assertEqual(len(table._cellvalues[0]),10)

    def test_control_only_does_not_run_checkm_workflow(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); genomes=root/'genomes'; genomes.mkdir()
            (genomes/'Control01.fasta').write_text('>lambda\nACGT\n')
            with patch('bactipipe.scripts.process_data.subprocess.run',return_value=SimpleNamespace(returncode=0,stdout=b'::: CheckM v1.2.4\n')) as run:
                result=checkM_stats(genomes,root/'checkm',logfile=root/'log',exclude_samples=['Control01'])
            self.assertEqual(result,({},'v1.2.4'))
            self.assertEqual(run.call_count,1)  # version query only
            self.assertTrue((genomes/'Control01.fasta').is_file())

    def test_checkm_mixed_batch_receives_only_bacterial_assemblies(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); genomes=root/'genomes'; genomes.mkdir()
            for sample in ('Control01','S1'):
                (genomes/f'{sample}.fasta').write_text('>seq\nACGT\n')
            def run_tool(args, **kwargs):
                if args[1]=='lineage_wf':
                    self.assertEqual([p.name for p in Path(args[-2]).glob('*.fasta')],['S1.fasta'])
                    self.assertTrue((Path(args[-2])/'S1.fasta').is_symlink())
                    return SimpleNamespace(returncode=1)  # no actual CheckM execution
                return SimpleNamespace(returncode=0,stdout=b'::: CheckM v1.2.4\n')
            with patch('bactipipe.scripts.process_data.subprocess.run',side_effect=run_tool), patch('bactipipe.scripts.process_data.time_print'):
                checkM_stats(genomes,root/'checkm',logfile=root/'log',exclude_samples=['Control01'])
            self.assertTrue((genomes/'Control01.fasta').is_file())

    def test_stale_bacterial_reconciliation_does_not_replace_control(self):
        from bactipipe.batch_species import update_summary_identifications
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'summary.tsv'
            path.write_text('Sample\tExpected organism\tIdentified organism\nControl01\tLambda\t'+CONTROL_IDENTIFICATION+'\n')
            original=path.read_text()
            update_summary_identifications(path,{'Control01':object()})
            self.assertEqual(path.read_text(),original)


if __name__=='__main__': unittest.main()
