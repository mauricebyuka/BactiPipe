from io import StringIO
from pathlib import Path
import tempfile
import unittest

from Bio import Phylo
from bactipipe.placement_visuals import prepare_tree, tip_label, ani_rows, render_tree, render_ani, shared_reference_panel


class PlacementVisualTests(unittest.TestCase):
    def test_either_af_plot_and_shared_panel_use_larger_fraction(self):
        evidence = {'identification_thresholds': {'minimum_ani': 95,
            'minimum_aligned_fraction': 50, 'aligned_fraction_mode': 'either'},
            'reference_audit': [dict(accession='GCF_1.1', species_name='Example species',
                ani=95.44, aligned_fraction_query=40, aligned_fraction_reference=71.42)]}
        rows, ani, af = ani_rows(evidence)
        self.assertEqual((ani, af, rows[0]['decision_af']), (95, 50, 71.42))
        panel = shared_reference_panel(evidence)
        self.assertEqual(panel[0]['decision_af'], 71.42)
        self.assertEqual(panel[0]['aligned_fraction_query'], 40)
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / 'ani.svg'
            render_ani(evidence, output, panel=panel)
            text = output.read_text()
            self.assertIn('Larger of query / reference aligned fractions (%)', text)
            self.assertIn('Highest maximum aligned fraction: R01 (71.42%)', text)

    def test_unnamed_species_hidden_only_in_shared_presentation(self):
        with tempfile.TemporaryDirectory() as tmp:
            tree=Path(tmp)/'tree.nwk'
            tree.write_text("(query:.1,'Named species [GCF_1.1]':.2,'Mannheimia sp041 [GCF_2.1]':.3);")
            unknown=dict(accession='GCF_2.1',species_name='Unresolved reference GCF_2.1',
                ani=83,aligned_fraction_query=3,aligned_fraction_reference=4)
            evidence={'reference_audit':[unknown]}
            panel=shared_reference_panel(evidence,tree)
            self.assertEqual([r['accession'] for r in panel],['GCF_1.1'])
            self.assertEqual(evidence['reference_audit'],[unknown])
            self.assertEqual(len(ani_rows(evidence)[0]),1)
    def test_shared_panel_joins_exact_accession_not_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            tree=Path(tmp)/'tree.nwk'
            tree.write_text("(query:.1,'Example species [GCF_1.1]':.2,'Other species [GCF_2.1]':.3);")
            evidence={'reference_audit':[
                dict(accession='GCF_1.2',species_name='Example species',ani=90,aligned_fraction_query=7,aligned_fraction_reference=8),
                dict(accession='GCF_2.1',species_name='Older species name',ani=89,aligned_fraction_query=5,aligned_fraction_reference=6)]}
            panel=shared_reference_panel(evidence,tree)
            by_acc={row['accession']:row for row in panel}
            self.assertEqual(by_acc['GCF_1.2']['availability'],'ANI comparison only')
            self.assertEqual(by_acc['GCF_1.1']['availability'],'ANI not measurable')
            self.assertEqual(by_acc['GCF_2.1']['species_name'],'Other species')
            self.assertEqual([row['panel_id'] for row in panel],['R01','R02','R03'])
            render_tree(tree,Path(tmp)/'tree.svg',panel=panel)
            render_ani(evidence,Path(tmp)/'ani.svg',panel=panel)
            tree_svg=(Path(tmp)/'tree.svg').read_text()
            ani_svg=(Path(tmp)/'ani.svg').read_text()
            self.assertIn('ANI comparison only',tree_svg)
            self.assertIn('Not measurable',ani_svg)
            self.assertIn('Highest minimum aligned fraction: R02 (7.00%)',ani_svg)
            self.assertIn('R03  Other species',tree_svg)
            self.assertNotIn('Older species name',ani_svg)

    def test_shared_panel_preserves_outgroup_and_bounded_references(self):
        with tempfile.TemporaryDirectory() as tmp:
            tree=Path(tmp)/'tree.nwk'
            tree.write_text('(query:.1,'+','.join(f"'Species {i} [GCF_{i}.1]':{i}" for i in range(1,15))+');')
            evidence={'reference_audit':[dict(accession=f'GCA_{i}.1',species_name=f'ANI {i}',
                ani=90,aligned_fraction_query=i,aligned_fraction_reference=i) for i in range(20,35)]}
            panel=shared_reference_panel(evidence,tree,outgroup={'reference_accession':'GCF_14.1'})
            self.assertEqual(len(panel),10)
            self.assertIn('GCF_14.1',[r['accession'] for r in panel])
            self.assertIn('GCF_1.1',[r['accession'] for r in panel])
            self.assertEqual(sum(r['availability']=='ANI comparison only' for r in panel),5)

    def test_pruning_preserves_distances_and_original(self):
        tree=Phylo.read(StringIO('(query:.1,'+','.join(f'R{i}:.{i+1}' for i in range(15))+');'),'newick')
        result,nearest=prepare_tree(tree,'query')
        self.assertEqual(len(result.get_terminals()),11)
        self.assertEqual(len(tree.get_terminals()),16)
        for tip in result.get_terminals():
            self.assertAlmostEqual(result.distance('query',tip.name),tree.distance('query',tip.name))
        self.assertEqual(nearest,'R0')

    def test_missing_query_and_invalid_lengths_rejected(self):
        for source in ('(not_query:.1,R:.2);','(query:-.1,R:.2);'):
            with self.assertRaises(ValueError): prepare_tree(Phylo.read(StringIO(source),'newick'),'query')

    def test_label_preserves_gtdb_suffix_separately(self):
        self.assertEqual(tip_label('Enterococcus cecorum (GTDB: Enterococcus_E cecorum_A) [GCA_023221775.1]'),
            ('Enterococcus cecorum','GCA_023221775.1','Enterococcus_E cecorum_A'))

    def test_both_afs_and_saved_thresholds_are_used(self):
        evidence={'identification_thresholds':{'minimum_ani':96,'minimum_aligned_fraction':70},
            'reference_audit':[dict(ani=97,aligned_fraction_query=80,aligned_fraction_reference=42,species_name='Example species')]}
        rows,ani,af=ani_rows(evidence)
        self.assertEqual((ani,af),(96,70))
        self.assertEqual(rows[0]['minimum_af'],42)

    def test_missing_measurements_are_not_zero(self):
        evidence={'reference_audit':[dict(ani=97,aligned_fraction_query=80)]}
        self.assertEqual(ani_rows(evidence)[0],[])

    def test_empty_reference_audit_preserves_original_measured_comparisons(self):
        row=dict(ani=84,aligned_fraction_query=3,aligned_fraction_reference=4)
        self.assertEqual(ani_rows({'reference_audit':[], 'closest_species':[row]})[0][0]['minimum_af'],3)

    def test_png_rendering(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'ani.png'
            render_ani({},path)
            self.assertTrue(path.read_bytes().startswith(b'\x89PNG\r\n\x1a\n'))

    def test_ten_reference_limit_and_distinct_cluster_preserved(self):
        evidence={'reference_audit':[dict(ani=97,aligned_fraction_query=80,aligned_fraction_reference=75,
            species_name='Enterococcus cecorum',gtdb_taxonomy=f's__Enterococcus_E cecorum_{i}',accession=f'GCA_{i}.1') for i in range(14)]}
        self.assertEqual(len(ani_rows(evidence)[0]),10)

    def test_render_does_not_invent_bootstrap_or_type_claims(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); tree=root/'tree.nwk'
            tree.write_text("((query:.01,'Example species [GCF_123.1]':.02)100:.01,'Different species [GCF_456.1]':.1);")
            original=tree.read_bytes()
            render_tree(tree,root/'tree.svg',sample_name='SAMPLE-A')
            svg=(root/'tree.svg').read_text()
            self.assertIn('SAMPLE-A',svg)
            self.assertNotIn('NEAREST DISPLAYED REFERENCE',svg)
            self.assertIn('Nearest reference',svg)
            self.assertNotIn('[TYPE]',svg)
            self.assertNotIn('>100<',svg)
            self.assertEqual(original,tree.read_bytes())
            render_ani({},root/'ani.svg',sample_name='SAMPLE-A')
            self.assertIn('No reportable ANI comparisons',(root/'ani.svg').read_text())


if __name__=='__main__': unittest.main()
