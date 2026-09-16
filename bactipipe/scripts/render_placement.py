"""Redraw retained tree/ANI evidence without repeating biological analysis."""
import argparse
import json
from pathlib import Path
from bactipipe.placement_visuals import render_tree, render_ani, shared_reference_panel


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tree',help='Saved compact Newick')
    p.add_argument('--query-name',default='query',help='Exact leaf identifier in Newick')
    p.add_argument('--sample',required=True,help='Displayed sample name')
    p.add_argument('--confirmation-json',help='Retained whole-assembly confirmation JSON')
    p.add_argument('--placement-json',help='Retained placement JSON (including outgroup metadata)')
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--format',choices=('svg','png','pdf'),default='svg')
    a=p.parse_args()
    if not a.tree and not a.confirmation_json: p.error('Supply --tree or --confirmation-json')
    a.output.mkdir(parents=True,exist_ok=True)
    evidence=json.loads(Path(a.confirmation_json).read_text()) if a.confirmation_json else {}
    placement=json.loads(Path(a.placement_json).read_text()) if a.placement_json else {}
    source_hash=placement.get('provenance',{}).get('assembly_sha256')
    ani_hash=evidence.get('provenance',{}).get('assembly_sha256')
    if source_hash and ani_hash and source_hash != ani_hash:
        p.error('Tree and ANI evidence belong to different assemblies')
    panel=shared_reference_panel(evidence,a.tree,query_name=a.query_name,outgroup=placement.get('outgroup'))
    if a.tree: render_tree(a.tree,a.output/f'tree.{a.format}',query_name=a.query_name,sample_name=a.sample,
                          panel=panel,outgroup=placement.get('outgroup'))
    if a.confirmation_json:
        render_ani(evidence,a.output/f'ani.{a.format}',sample_name=a.sample,panel=panel)


if __name__=='__main__': main()
