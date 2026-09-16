from __future__ import annotations

import tempfile
import unittest
import importlib.util
import gzip
from pathlib import Path
from unittest.mock import patch

from bactipipe.phylogenetic_placement import (
    _compact_tree,
    _find_query_tree,
    _load_species_names,
    _distances_from_query,
    _induced_tree,
    render_bac120_placement,
)


class PhylogeneticPlacementTests(unittest.TestCase):
    def test_ani_nominee_uses_exact_accession_from_full_tree(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); source=root/'full.nwk'; metadata=root/'metadata.tsv'
            source.write_text('(query:1,RS_GCF_1.1:1,RS_GCF_2.1:2,RS_GCF_3.1:3,RS_GCF_4.1:4);')
            metadata.write_text('accession\tgtdb_taxonomy\n'+''.join(
                f'RS_GCF_{i}.1\ts__Species {i}\n' for i in range(1,5)))
            with patch('bactipipe.phylogenetic_placement._render_svg'):
                result=_compact_tree(source_tree=source,query_name='query',metadata_path=metadata,
                    output_directory=root,maximum_species=3,preferred_accessions=['GCF_4.1','GCF_2.2'])
            saved=(root/result.tree_newick).read_text()
            self.assertIn('GCF_1.1',saved)
            self.assertIn('GCF_4.1',saved)
            self.assertNotIn('GCF_2.2',saved)

    def test_outgroup_is_outside_full_tree_ingroup_and_distances_preserved(self):
        from Bio import Phylo
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            source=root/'full.nwk'
            source.write_text('(((query:.01,RS_GCF_1.1:.02):.01,RS_GCF_2.1:.03):.1,RS_GCF_3.1:.2);')
            metadata=root/'metadata.tsv'
            metadata.write_text('accession\tgtdb_taxonomy\n'+''.join(
                f'RS_GCF_{i}.1\ts__Species {i}\n' for i in range(1,4)))
            with patch('bactipipe.phylogenetic_placement._render_svg'):
                result=_compact_tree(source_tree=source,query_name='query',metadata_path=metadata,
                    output_directory=root,maximum_species=3)
            self.assertEqual(result.outgroup['reference_accession'],'GCF_3.1')
            compact=Phylo.read(str(root/result.tree_newick),'newick')
            original=Phylo.read(str(source),'newick')
            self.assertAlmostEqual(compact.distance('query','Species 3 [GCF_3.1]'),
                                   original.distance('query','RS_GCF_3.1'))
            self.assertEqual(len(compact.get_terminals()),4)

    def test_root_polytomy_does_not_invent_outgroup(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); source=root/'full.nwk'; metadata=root/'metadata.tsv'
            source.write_text('(query:1,RS_GCF_1.1:1,RS_GCF_2.1:1,RS_GCF_3.1:2);')
            metadata.write_text('accession\tgtdb_taxonomy\n'+''.join(
                f'RS_GCF_{i}.1\ts__Species {i}\n' for i in range(1,4)))
            with patch('bactipipe.phylogenetic_placement._render_svg'):
                result=_compact_tree(source_tree=source,query_name='query',metadata_path=metadata,
                    output_directory=root,maximum_species=3)
            self.assertEqual(result.outgroup,{})

    def test_query_must_match_a_leaf_not_an_accession_substring(self):
        with tempfile.TemporaryDirectory() as directory:
            tree = Path(directory) / "placement.bac120.classify.tree"
            tree.write_text("(query_other:0.1,reference_query:0.2);\n")
            self.assertIsNone(_find_query_tree(Path(directory), "query"))

    def test_species_suffix_is_not_silently_renamed(self):
        with tempfile.TemporaryDirectory() as directory:
            metadata = Path(directory) / "metadata.tsv"
            metadata.write_text("accession\tgtdb_taxonomy\nRS_GCF_1.1\td__Bacteria;s__Enterococcus_B cecorum\n")
            self.assertEqual(_load_species_names(metadata, {"GCF_1.1"})["GCF_1.1"], "Enterococcus_B cecorum")

    def test_distinct_clusters_with_same_ncbi_name_remain_distinct(self):
        with tempfile.TemporaryDirectory() as directory:
            metadata = Path(directory) / 'names.tsv'
            metadata.write_text('accession\tncbi_species\tgtdb_taxonomy\n'
                'RS_GCF_1.1\tEnterococcus cecorum\td__Bacteria;s__Enterococcus_E cecorum\n'
                'GB_GCA_2.1\tEnterococcus cecorum\td__Bacteria;s__Enterococcus_E cecorum_A\n')
            names = _load_species_names(metadata, {'GCF_1.1', 'GCA_2.1'})
            self.assertEqual(len(set(names.values())), 2)
            self.assertTrue(all('Enterococcus cecorum' in value for value in names.values()))

    def test_induced_tree_preserves_distances_and_source(self):
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO("(((query:1,A:2):3,B:4):5,(C:6,D:7):8);"), "newick")
        query = next(tip for tip in tree.get_terminals() if tip.name == "query")
        distances = _distances_from_query(tree, query)
        for tip in tree.get_terminals():
            self.assertAlmostEqual(distances[tip], tree.distance(query, tip))
        compact = _induced_tree(tree, {"query", "B", "D"})
        for first, second in (("query", "B"), ("query", "D"), ("B", "D")):
            self.assertAlmostEqual(tree.distance(first, second), compact.distance(first, second))
        self.assertEqual(len(tree.get_terminals()), 5)

    def test_existing_batch_can_be_rendered_without_rerunning_gtdb(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "placement.bac120.classify.tree").write_text("(query:1,RS_GCF_1.1:2,RS_GCF_2.1:3);")
            metadata = root / "metadata.tsv"
            metadata.write_text("accession\tgtdb_taxonomy\nRS_GCF_1.1\ts__Species alpha\nRS_GCF_2.1\ts__Species beta\n")
            (root / "placement.bac120.markers_summary.tsv").write_text(
                "name\tnumber_unique_genes\tnumber_multiple_genes\tnumber_missing_genes\nquery\t110\t2\t8\n")
            with patch("bactipipe.phylogenetic_placement._render_svg") as render, \
                    patch("bactipipe.phylogenetic_placement.subprocess.run") as run:
                for _ in range(2):
                    result = render_bac120_placement(gtdb_output=root, query_name="query",
                        metadata_path=metadata, output_directory=root / "compact")
                    self.assertEqual(result.status, "available")
                    self.assertEqual(result.marker_evidence["single_copy_markers"], 110)
                    self.assertEqual(result.marker_evidence["multiple_copy_markers"], 2)
                    self.assertIn("not a confirmed", result.interpretation)
                run.assert_not_called()
                self.assertEqual(render.call_count, 2)

    def test_no_tree_retains_marker_evidence_not_a_synthetic_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = render_bac120_placement(gtdb_output=root, query_name="missing",
                metadata_path=root / "absent.tsv", output_directory=root / "result")
            self.assertEqual(result.status, "unavailable")
            self.assertNotIn("single_copy_markers", result.marker_evidence)

    def test_gtdb_metadata_can_remain_compressed(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            metadata = Path(temporary_dir) / "metadata.tsv.gz"
            with gzip.open(metadata, "wt", newline="") as handle:
                handle.write(
                    "accession\tncbi_species\tgtdb_taxonomy\n"
                    "RS_GCF_1.1\ts__Species alpha\td__Bacteria;s__Species_alpha\n"
                )
            self.assertEqual(
                _load_species_names(metadata, {"GCF_1.1"}),
                {"GCF_1.1": "Species alpha (GTDB: Species_alpha)"},
            )

    def test_current_gtdb_ncbi_taxonomy_column_supplies_species_name(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            metadata = Path(temporary_dir) / "metadata.tsv"
            metadata.write_text(
                "accession\tncbi_taxonomy\tgtdb_taxonomy\n"
                "RS_GCF_1.1\td__Bacteria;p__P;c__C;o__O;f__F;g__G;"
                "s__Veterinary species\td__Bacteria;s__GTDB_species\n"
            )
            self.assertEqual(
                _load_species_names(metadata, {"GCF_1.1"}),
                {"GCF_1.1": "Veterinary species (GTDB: GTDB_species)"},
            )

    def test_finds_standard_unsuffixed_classification_tree(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            tree = Path(temporary_dir) / "placement.bac120.classify.tree"
            tree.write_text("(query:0.1,reference:0.2);\n")
            self.assertEqual(_find_query_tree(Path(temporary_dir), "query"), tree)

    @unittest.skipUnless(importlib.util.find_spec("Bio"), "Biopython is not installed")
    def test_compact_tree_keeps_one_reference_per_species(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            tree = root / "placement.bac120.classify.tree.1.tree"
            tree.write_text(
                "((query:0.01,RS_GCF_1.1:0.02):0.01,"
                "(RS_GCF_2.1:0.03,RS_GCF_3.1:0.04):0.02);\n"
            )
            metadata = root / "metadata.tsv"
            metadata.write_text(
                "accession\tncbi_species\tgtdb_taxonomy\n"
                "RS_GCF_1.1\ts__Species alpha\td__Bacteria;s__Species_alpha\n"
                "RS_GCF_2.1\ts__Species beta\td__Bacteria;s__Species_beta\n"
                "RS_GCF_3.1\ts__Species beta\td__Bacteria;s__Species_beta\n"
            )

            def fake_render(_tree, destination, **_kwargs):
                destination.write_text("<svg/>\n")

            with patch(
                "bactipipe.phylogenetic_placement._render_svg",
                side_effect=fake_render,
            ):
                result = _compact_tree(
                    source_tree=tree,
                    query_name="query",
                    metadata_path=metadata,
                    output_directory=root,
                    maximum_species=10,
                )
            self.assertEqual(result.status, "available")
            self.assertEqual(result.closest_species, ("Species alpha (GTDB: Species_alpha)", "Species beta (GTDB: Species_beta)"))
            neighbors = (root / "closest_species_tree_neighbors.tsv").read_text()
            self.assertEqual(neighbors.count("Species beta"), 1)
            self.assertTrue((root / "closest_species_tree.nwk").exists())


if __name__ == "__main__":
    unittest.main()
