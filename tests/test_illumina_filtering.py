import json
import random
import shutil
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from bactipipe.illumina_filtering import (
    select_reads, needs_filtering, combined_metrics, confident_depth_organism,
    refresh_depth_row, assemble_selected,
)
from bactipipe.scripts.qc_plan import illumina_minimum_coverage


def write_reads(path, qualities):
    rng = random.Random(12)
    with path.open('w') as out:
        for index, quality in enumerate(qualities):
            seq = ''.join(rng.choice('ACGT') for _ in quality)
            out.write(f'@read{index}\n{seq}\n+\n{quality}\n')


class IlluminaFilteringTests(unittest.TestCase):
    def test_trigger_boundaries_and_depth_tiers(self):
        self.assertTrue(needs_filtering(28.999999, 28))
        self.assertFalse(needs_filtering(29, 28))
        self.assertFalse(needs_filtering(30, 28))
        self.assertTrue(needs_filtering(30.9, 30))
        self.assertFalse(needs_filtering(31, 30))
        for q, minimum in ((28,50),(28.99999,50),(29,40),(29.996,40),(30,30)):
            self.assertEqual(illumina_minimum_coverage(q,None),minimum)
            self.assertEqual(illumina_minimum_coverage(q,'Escherichia coli'),minimum+10)
        with self.assertRaises(ValueError): needs_filtering(float('nan'),28)

    def test_combined_mean_is_base_weighted(self):
        q=SimpleNamespace(combined_total_bases=1000,combined_num_reads=2,combined_quality_sum=100*40+900*27)
        metrics=combined_metrics(q)
        self.assertEqual(metrics['mean_quality'],28.3)
        self.assertTrue(needs_filtering(metrics['mean_quality'],28))

    def test_only_confident_species_can_add_ecoli_requirement(self):
        def fields(match,name='Escherichia coli'):
            return SimpleNamespace(match_percent=match,coverage_organism=name,identified_organism=name)
        self.assertIsNone(confident_depth_organism(fields('N/A')))
        self.assertEqual(confident_depth_organism(fields('99.8')),'Escherichia coli')
        self.assertEqual(confident_depth_organism(fields('99.8','Salmonella enterica')),'Salmonella enterica')

    def test_raw_q29_skips_fastp_and_preserves_full_precision_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); reads=root/'reads.fastq'; qc=root/'qc'; qc.mkdir()
            # 29.996 rounds to 30.00 for display, but must keep the 40X tier.
            write_reads(reads, ['?'*996+'>'*4])
            with patch('bactipipe.illumina_filtering.fastp_version',side_effect=AssertionError('fastp must not run')):
                result=select_reads([str(reads)],qc,'S',28,1,str(root/'log'))
            self.assertFalse(result['filter_applied'])
            self.assertEqual(result['selected_paths'],[str(reads)])
            self.assertEqual(result['selected']['mean_quality'],29.996)
            with patch('bactipipe.scripts.qualityProc.ProcQuality',side_effect=AssertionError('cached metrics should be reused')):
                reused=select_reads([str(reads)],qc,'S',28,1,str(root/'log'))
            self.assertEqual(reused['selected']['mean_quality'],29.996)
            self.assertEqual(illumina_minimum_coverage(reused['selected']['mean_quality'],None),40)

    @unittest.skipUnless(shutil.which('fastp'), 'fastp required')
    def test_real_filtering_uses_both_mates_and_selected_yield(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); qc=root/'qc'; qc.mkdir()
            paths=[root/'R1.fastq',root/'R2.fastq']
            for path in paths:
                write_reads(path,[':'*150]*80+['I'*150]*20) # Q25 and Q40
            result=select_reads([str(p) for p in paths],qc,'S',28,2,str(root/'log'))
            self.assertTrue(result['filter_applied'])
            self.assertEqual(result['raw']['mean_quality'],28)
            self.assertEqual(result['selected']['mean_quality'],40)
            self.assertLess(result['selected']['total_bases'],result['raw']['total_bases'])
            self.assertEqual(len(result['selected_paths']),2)
            self.assertNotEqual(result['selected_paths'],[str(p) for p in paths])
            self.assertIn('--disable_length_filtering',result['fastp_command'])
            for p in paths: self.assertTrue(p.exists())

    @unittest.skipUnless(shutil.which('fastp'), 'fastp required')
    def test_no_reads_after_filtering_does_not_fall_back_to_raw(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); qc=root/'qc'; qc.mkdir(); reads=root/'reads.fastq'
            write_reads(reads,['&'*150]*10)
            with self.assertRaisesRegex(ValueError,'No usable reads'):
                select_reads([str(reads)],qc,'S',28,1,str(root/'log'))
            self.assertFalse((qc/'S'/'illumina_read_selection.json').exists())

    def test_malformed_fastq_does_not_publish_partial_metrics(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); qc=root/'qc'; qc.mkdir(); reads=root/'reads.fastq'
            reads.write_text('@r\nACGT\n+\nIII\n')
            with self.assertRaisesRegex(ValueError,'Malformed'):
                select_reads([str(reads)],qc,'S',28,1,str(root/'log'))

    def test_reconciliation_uses_full_precision_and_assembly_not_expectation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); audit=root/'qc_out'/'S'; audit.mkdir(parents=True)
            (audit/'illumina_read_selection.json').write_text(json.dumps({'selected':{'mean_quality':29.996,'total_bases':18000}}))
            genomes=root/'assemblies'/'genomes'; genomes.mkdir(parents=True)
            (genomes/'S.fasta').write_text('>a\n'+'A'*400+'\n')
            row={'Sample':'S','Expected organism':'Salmonella enterica','min_cov':'30','cov_verdict':'Pass','Coverage Depth':'999X'}
            confident=SimpleNamespace(match_percent='99',coverage_organism='Escherichia coli',identified_organism='Escherichia coli')
            refresh_depth_row(row,confident,root/'qc_out'/'temp_qc_summary.tsv')
            self.assertEqual((row['min_cov'],row['Coverage Depth'],row['cov_verdict']),('50','45.00X','Fail'))
            confident.match_percent='N/A'
            refresh_depth_row(row,confident,root/'qc_out'/'temp_qc_summary.tsv')
            self.assertEqual((row['min_cov'],row['cov_verdict']),('40','Pass'))

    def test_assembly_cache_requires_selection_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            selection={'selection_id':'first','selected':{'total_bases':1000},'selected_paths':['selected.fastq']}
            def fake_assemble(**kwargs):
                out=Path(kwargs['assembly_dir'])/'genomes'; out.mkdir(parents=True)
                (out/'S.fasta').write_text('>new\nAAAA\n')
                self.assertEqual(kwargs['reads'],['selected.fastq'])
            with patch('bactipipe.scripts.process_data.assemble',side_effect=fake_assemble) as assemble:
                assemble_selected('S',selection,root,'Skesa','v1',1,None)
                assemble_selected('S',selection,root,'Skesa','v1',1,None)
                self.assertEqual(assemble.call_count,1)
                selection['selection_id']='changed'
                assemble_selected('S',selection,root,'Skesa','v1',1,None)
                self.assertEqual(assemble.call_count,2)

    def test_selective_retry_replaces_plots_and_assembly_identity(self):
        from bactipipe.scripts.sample_update import _replace_sample_outputs
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fresh, merged = root / 'fresh', root / 'merged'
            relatives = ('qc_out/S_metrics.png', 'qc_out/S_quality_distribution.png',
                         'assemblies/genomes/S.input.json')
            for relative in relatives:
                for folder, value in ((fresh, 'new'), (merged, 'old')):
                    path = folder / relative
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text(value)
            untouched = merged / 'qc_out' / 'OTHER_metrics.png'
            untouched.write_text('unrelated')
            _replace_sample_outputs(fresh, merged, 'S')
            for relative in relatives:
                self.assertEqual((merged / relative).read_text(), 'new')
            self.assertEqual(untouched.read_text(), 'unrelated')

    def test_runner_has_no_catalog_depth_gate(self):
        import ast
        source=(Path(__file__).parents[1]/'bactipipe/scripts/run_illumina_qc.py').read_text()
        self.assertNotIn('coverage_and_verdict(',source)
        self.assertNotIn('if resolution.is_unknown and identified_org',source)
        tree=ast.parse(source)
        calls=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name)]
        self.assertTrue(any(n.func.id=='prepare_batch' for n in calls))
        worker=(Path(__file__).parents[1]/'bactipipe/illumina_parallel.py').read_text()
        calls=[n for n in ast.walk(ast.parse(worker)) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name)]
        self.assertTrue(any(n.func.id=='assemble_selected' for n in calls))
        self.assertTrue(any(n.func.id=='select_reads' for n in calls))


if __name__=='__main__': unittest.main()
