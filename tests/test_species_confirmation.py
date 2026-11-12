import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from bactipipe.species_confirmation import evaluate_ani, investigate_assembly, reference_name, screen_assembly
from bactipipe.confirmation_summary import update_row


def hit(acc='GCA_023221775.1', ani=97.07, query=71.77, ref=79.12):
    return dict(Ref_file=acc+'.fna', ANI=ani, Align_fraction_query=query, Align_fraction_ref=ref)


def label(name='Enterococcus cecorum', taxid='533'):
    return dict(ncbi_taxonomy='d__Bacteria;s__'+name, ncbi_species_taxid=taxid,
        ncbi_organism_name=name, gtdb_taxonomy='s__Enterococcus_E cecorum_A', ncbi_genome_category='none')


class ConfirmationTests(unittest.TestCase):
    def evaluate(self, rows=None, labels=None, **kwargs):
        return evaluate_ani(rows or [hit()], labels or {'GCA_023221775.1': label()}, **kwargs)

    def test_cecorum_uses_ncbi_name_retains_cluster_evidence(self):
        result = self.evaluate()
        self.assertEqual(result['status'], 'confirmed_species')
        self.assertEqual(result['reported_organism'], 'Enterococcus cecorum')
        self.assertEqual(result['match_percent'], 97.07)
        self.assertIn('cecorum_A', result['reference_audit'][0]['gtdb_taxonomy'])

    def test_either_aligned_fraction_required_without_rounding(self):
        for row in (hit(query=49.99, ref=49.99), hit(ani=94.99)):
            self.assertEqual(self.evaluate([row])['status'], 'no_confident_identification')
        for row in (hit(query=50, ref=40), hit(query=40, ref=50), hit(ani=95)):
            result = self.evaluate([row])
            self.assertEqual(result['status'], 'confirmed_species')
            self.assertEqual(result['identification_thresholds']['aligned_fraction_mode'], 'either')

    def test_historical_both_af_policy_is_not_reinterpreted(self):
        evidence = {'identification_thresholds': {'minimum_ani': 95, 'minimum_aligned_fraction': 65}}
        self.assertEqual(self.evaluate([hit(query=62.49, ref=71.42)], evidence=evidence)['status'],
                         'no_confident_identification')
        evidence['identification_thresholds'].update(minimum_aligned_fraction=50, aligned_fraction_mode='either')
        self.assertEqual(self.evaluate([hit(query=62.49, ref=71.42)], evidence=evidence)['status'],
                         'confirmed_species')

    def test_low_alignment_evidence_is_retained_but_never_identifies(self):
        result=self.evaluate([hit(ani=99,query=2,ref=3)])
        self.assertEqual(result['status'],'no_confident_identification')
        self.assertIsNone(result['match_percent'])
        self.assertEqual(result['reference_audit'][0]['aligned_fraction_query'],2)

    def test_thresholds_from_original_analysis_preserved(self):
        result = self.evaluate(evidence={'identification_thresholds': {'minimum_ani': 98}})
        self.assertEqual(result['status'], 'no_confident_identification')

    def test_competing_species_remains_ambiguous(self):
        result = self.evaluate([hit(), hit('GCA_123456789.1', ani=97.0)],
            {'GCA_023221775.1': label(), 'GCA_123456789.1': label('Enterococcus faecium', '1352')})
        self.assertEqual(result['status'], 'ambiguous_species')

    def test_prior_competitor_not_erased(self):
        competitor = dict(accession='old', species_taxid='1352', species_name='Enterococcus faecium',
            reference_path='old.fasta', ani=97.0, aligned_fraction_reference=79, aligned_fraction_query=72)
        self.assertEqual(self.evaluate(evidence={'closest_species': [competitor]})['status'], 'ambiguous_species')

    def test_noncallable_reference_never_confirms(self):
        for metadata in ({}, label('Enterococcus sp.'), dict(label(), ncbi_species_taxid='')):
            result = self.evaluate(labels={'GCA_023221775.1': metadata})
            self.assertEqual(result['status'], 'no_confident_identification')

    def test_named_gtdb_mag_and_environmental_references_can_confirm(self):
        for category in ('derived from metagenome', 'derived from environmental sample', 'none'):
            with self.subTest(category=category):
                metadata = dict(label(), ncbi_genome_category=category)
                result = self.evaluate(labels={'GCA_023221775.1': metadata})
                self.assertEqual(result['status'], 'confirmed_species')
                self.assertEqual(result['reported_organism'], 'Enterococcus cecorum')
                self.assertEqual(result['reference_audit'][0]['ncbi_genome_category'], category)
                self.assertTrue(result['reference_audit'][0]['callable'])

    def test_p2629621_mag_and_isolate_agree_on_species(self):
        result = self.evaluate([
            hit('GCA_946499615.1', ani=96.87, query=66.29, ref=76.76),
            hit('GCF_052119425.1', ani=95.44, query=62.49, ref=71.42)], {
            'GCA_946499615.1': dict(label(), ncbi_genome_category='derived from metagenome'),
            'GCF_052119425.1': label()})
        self.assertEqual(result['status'], 'confirmed_species')
        self.assertEqual(result['reported_organism'], 'Enterococcus cecorum')
        self.assertEqual(result['match_percent'], 96.87)

    def test_mag_origin_does_not_bypass_numeric_or_competitor_gates(self):
        metadata = dict(label(), ncbi_genome_category='derived from metagenome')
        for row in (hit(ani=94.99), hit(query=49.99, ref=49.99)):
            self.assertEqual(self.evaluate([row], {'GCA_023221775.1': metadata})['status'], 'no_confident_identification')
        result = self.evaluate([hit(), hit('GCA_123456789.1', ani=97.0)], {
            'GCA_023221775.1': metadata, 'GCA_123456789.1': label('Enterococcus faecium', '1352')})
        self.assertEqual(result['status'], 'ambiguous_species')

    def test_generic_equi_not_promoted_to_subspecies(self):
        self.assertFalse(reference_name(label('Streptococcus equi', '1336'), 'x')[2])
        row = dict(label('Streptococcus equi', '1336'), ncbi_organism_name='Streptococcus equi subsp. zooepidemicus strain ABC')
        self.assertEqual(reference_name(row, 'x')[:2], ('Streptococcus equi subsp. zooepidemicus', '40041'))

    def test_tree_skipped_for_confirmed_ani(self):
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, 'confirmation.json').write_text(json.dumps(self.evaluate()))
            with patch('bactipipe.species_confirmation.screen_assembly', return_value={'directory':tmp, 'status':'identified'}), patch('bactipipe.species_confirmation.place_iteratively') as tree:
                result = investigate_assembly(assembly='query')
                self.assertEqual(result['status'], 'identified')
                tree.assert_not_called()

    def test_tree_only_after_successful_unresolved_screen(self):
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, 'confirmation.json').write_text(json.dumps(self.evaluate([hit(ani=93)])))
            Path(tmp, 'phylogenetic_placement.json').write_text('{"status":"available"}')
            with patch('bactipipe.species_confirmation.screen_assembly', return_value={'directory':tmp, 'status':'unresolved'}), patch('bactipipe.species_confirmation.place_iteratively', return_value={'directory':tmp, 'status':'available'}) as tree:
                self.assertEqual(investigate_assembly(assembly='query')['status'], 'available')
                tree.assert_called_once()

    def test_technical_failure_does_not_start_tree(self):
        with patch('bactipipe.species_confirmation.screen_assembly', side_effect=RuntimeError('failed')), patch('bactipipe.species_confirmation.place_iteratively') as tree:
            with self.assertRaises(RuntimeError): investigate_assembly(assembly='query')
            tree.assert_not_called()

    def test_invalid_ani_rejected(self):
        with self.assertRaises(ValueError): self.evaluate([hit(ani=float('nan'))])

    def test_confirmation_does_not_override_failed_qc(self):
        row = {'Expected organism': 'Enterococcus cecorum', 'Mean quality': '25',
            'Coverage Depth': '100X', 'Required depth':'50X', 'CheckM completeness':'99',
            'CheckM contamination':'0', 'Overall Quality':'Review'}
        self.assertEqual(update_row(row, self.evaluate())['Overall Quality'], 'Pass')
        self.assertEqual(update_row(dict(row, **{'Mean quality': '9'}), self.evaluate())['Overall Quality'], 'Fail')
        self.assertEqual(update_row(dict(row, **{'CheckM completeness': 'N/A'}), self.evaluate())['Overall Quality'], 'Review')
        self.assertEqual(update_row(dict(row, **{'Expected organism':'Salmonella enterica'}), self.evaluate())['Overall Quality'], 'Pass')

    def test_new_unresolved_confirmation_clears_legacy_match_preserves_measurements(self):
        row = {'Expected organism': 'Unknown', 'Identified organism': 'Staphylococcus aureus',
            '% Match': '96.40', 'Mean quality': '25', 'Coverage Depth': '100X',
            'Required depth': '50X', 'CheckM completeness': '99', 'CheckM contamination': '0'}
        result = {'status': 'no_confident_identification'}
        updated = update_row(row, result)
        self.assertEqual(updated['Identified organism'], 'No confident identification')
        self.assertEqual(updated['% Match'], 'N/A')
        self.assertEqual(updated['Coverage Depth'], '100X')
        self.assertEqual(updated['Overall Quality'], 'Review')
        self.assertEqual(row['% Match'], '96.40')
        self.assertEqual(update_row(dict(row, **{'CheckM contamination': '79.51'}), result)['Overall Quality'], 'Fail')

    def test_ecoli_confirmation_recalculates_illumina_depth_requirement(self):
        result = dict(self.evaluate(), reported_organism='Escherichia coli')
        row = {'Expected organism': 'Escherichia coli', 'Mean quality':'30', 'Coverage Depth':'35X',
            'Required depth':'30X', 'CheckM completeness':'99', 'CheckM contamination':'0'}
        updated = update_row(row, result, platform='illumina', quality=28, selected_mean_quality=30.01)
        self.assertEqual(updated['Required depth'], '40X')
        self.assertEqual(updated['Coverage Depth'], '35X')
        self.assertEqual(updated['Overall Quality'], 'Fail')

    def test_ani_reuse_requires_identical_entire_assembly(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            database = root / 'data/skani/database'; database.mkdir(parents=True)
            for name in ('metadata.tsv','index.db','sketches.db','markers.bin'):
                (database / name).write_text('immutable reference')
            (root / 'MD5SUM.txt').write_text('verified package')
            (root / 'data/metadata').mkdir()
            (root / 'data/metadata/metadata.txt').write_text('VERSION_DATA=r232')
            metadata = root / 'labels.tsv'; metadata.write_text('labels')
            assembly = root / 'query.fasta'; assembly.write_text('>contig\nAAAA\n')
            def run(command, *args):
                self.assertEqual(command[command.index('--min-af')+1],'1')
                destination = Path(command[command.index('-o')+1])
                destination.write_text('Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\n'
                    f'GCA_023221775.1.fna\t{assembly}\t97.07\t79.12\t71.77\n')
            options = dict(assembly=assembly, reference_root=root/'data', cache_root=root/'cache', metadata=metadata)
            with patch('bactipipe.species_confirmation.subprocess.check_output', return_value='skani 0.3.2'), patch('bactipipe.species_confirmation.reference_labels', return_value={'GCA_023221775.1':label()}), patch('bactipipe.species_confirmation.execute', side_effect=run) as engine:
                first = screen_assembly(**options)
                self.assertEqual(first['status'], 'identified')
                self.assertEqual(screen_assembly(**options)['reuse'], 'assembly')
                self.assertEqual(engine.call_count, 1)
                # Even a FASTA-only change must repeat ANI; marker equality is
                # irrelevant to the whole-assembly confirmation cache.
                assembly.write_text('>different_header\nAAAA\n')
                second = screen_assembly(**options)
                self.assertNotEqual(first['directory'], second['directory'])
                self.assertEqual(engine.call_count, 2)


if __name__ == '__main__':
    unittest.main()
