import json
from pathlib import Path

import pytest

from bactipipe.mlst_catalog import resolve_mlst, species_rows
from bactipipe.scripts import install_mlst_catalog as installer
from bactipipe.scripts import type_genomes as typing


def installed(root, scheme='bsubtilis_test', species='Bacillus subtilis'):
    entry = dict(scheme=scheme, species=[species], source='https://rest.pubmlst.org/db/test/schemes/1',
                 status='installed', loci=['a', 'b'])
    folder = root / scheme
    folder.mkdir(parents=True)
    (folder / (scheme+'.fsa')).write_text('>a_1\nACGT\n>b_1\nACGT\n')
    (folder / (scheme+'.tsv')).write_text('ST\ta\tb\n1\t1\t1\n')
    entry['files'] = installer.validate_scheme(root, entry)
    (root / 'mlst_catalog.json').write_text(json.dumps(dict(schema_version=1,
        schemes=[entry], species=species_rows([entry]))))
    return entry


def test_exact_species_lookup_and_integrity(tmp_path):
    installed(tmp_path)
    assert resolve_mlst('  BACILLUS   subtilis ', tmp_path)['scheme'] == 'bsubtilis_test'
    for name in ('Bacillus', 'Bacillus subtiliss', 'Bacillus subtilis subsp. spizizenii'):
        assert resolve_mlst(name, tmp_path)['status'] == 'no_matching_scheme'
    (tmp_path/'bsubtilis_test/bsubtilis_test.tsv').write_text('changed')
    with pytest.raises(ValueError, match='has changed'):
        resolve_mlst('Bacillus subtilis', tmp_path)


def test_configured_workflows_always_win(tmp_path):
    (tmp_path/'mlst_catalog.json').write_text('invalid')
    assert resolve_mlst('ecoli', tmp_path)['scheme'] == 'ecoli'
    assert resolve_mlst('klebsiella', tmp_path)['scheme'] == 'kpneumoniae'
    assert resolve_mlst('salmonella', tmp_path)['scheme'] == 'senterica'


def test_multiple_schemes_require_explicit_preference():
    a = dict(scheme='a', species=['Example species'], source='url1')
    b = dict(scheme='b', species=['Example species'], source='url2')
    assert species_rows([a,b])[0]['preferred_scheme'] is None
    rows = species_rows([a,b], [dict(species='Example species',source='url2')])
    assert rows[0]['preferred_scheme'] == 'b'
    # Never silently switch to another scheme if the preferred one disappears.
    assert species_rows([a], [dict(species='Example species',source='url2')])[0]['preferred_scheme'] is None


def test_missing_catalog_is_explicit(tmp_path):
    assert resolve_mlst('Streptococcus agalactiae',tmp_path)['status'] == 'catalog_unavailable'


def test_unavailable_preferred_scheme_is_not_substituted(tmp_path):
    entry=installed(tmp_path)
    entry['status']='unavailable'
    (tmp_path/'mlst_catalog.json').write_text(json.dumps(dict(schema_version=1,
        schemes=[entry],species=species_rows([entry]))))
    result=resolve_mlst('Bacillus subtilis',tmp_path)
    assert result['status']=='scheme_unavailable'
    assert result['scheme'] is None


def test_discovery_supplements_xml_and_excludes_whole_genome(monkeypatch):
    base='https://rest.pubmlst.org/db/test_seqdef'
    responses={installer.XML_URL:b'<data><species>Example species<mlst><database><profiles><url>https://rest.pubmlst.org/db/test_seqdef/schemes/1/profiles_csv</url></profiles></database></mlst></species></data>',
        installer.API_URL:json.dumps([dict(databases=[dict(name='test_seqdef',href=base)])]).encode(),
        base+'/schemes':json.dumps(dict(schemes=[
            dict(description='MLST',scheme=base+'/schemes/2'),
            dict(description='cgMLST',scheme=base+'/schemes/3'),
            dict(description='wgMLST',scheme=base+'/schemes/4')])).encode()}
    monkeypatch.setattr(installer,'MLST_SCHEMES',())
    rows,errors=installer.discover(responses.__getitem__)
    assert not errors
    assert {r['source'] for r in rows}=={base+'/schemes/1',base+'/schemes/2'}


def test_metadata_maps_exact_species_not_genus_or_complex():
    metadata = dict(loci=['https://rest.pubmlst.org/db/x/loci/a']*7,
                    profiles_csv='https://rest.pubmlst.org/db/x/schemes/1/profiles_csv',
                    description='MLST',ncbi_taxon=[dict(rank='genus',name='Bacillus')])
    row=installer.describe(dict(source='https://rest.pubmlst.org/db/x_seqdef/schemes/1',
        labels=['Bacillus subtilis#1','Bacillus spp.','Bacillus cereus complex']), lambda url:json.dumps(metadata).encode())
    assert row['species']==['Bacillus subtilis']
    metadata['description']='cgMLST'
    with pytest.raises(ValueError,match='conventional'):
        installer.describe(dict(source=row['source'],labels=[]),lambda url:json.dumps(metadata).encode())


def test_build_installs_and_records_unavailable_without_touching_seed(tmp_path, monkeypatch):
    source='https://rest.pubmlst.org/db/pubmlst_bsubtilis_seqdef/schemes/1'
    urls=[f'https://rest.pubmlst.org/db/pubmlst_bsubtilis_seqdef/loci/{n}' for n in ('a','b')]
    responses={source:json.dumps(dict(description='MLST',loci=urls,profiles_csv=source+'/profiles_csv',
        message='Public snapshot restriction',ncbi_taxon=[dict(rank='species',name='Bacillus subtilis')])).encode(),
        source+'/profiles_csv':b'ST\ta\tb\n1\t1\t1\n',
        urls[0]+'/alleles_fasta':b'>a_1\nACGT\n',urls[1]+'/alleles_fasta':b'>b_1\nACGT\n'}
    monkeypatch.setattr(installer,'discover',lambda fetch:([dict(source=source,labels=['Bacillus subtilis'])],[]))
    monkeypatch.setattr(installer.Downloads,'__call__',lambda self,url:responses[url])
    output,work=tmp_path/'version',tmp_path/'build'
    status=installer.build(output,work)
    assert status['installed']==1
    decision=resolve_mlst('Bacillus subtilis',output)
    assert decision['status']=='selected'
    assert decision['access_note']=='Public snapshot restriction'
    assert (output/'config').read_text().count(decision['scheme'])==1
    with pytest.raises(ValueError,match='exists'):
        installer.build(output,work)


def test_output_cannot_overlap_workspace_or_seed(tmp_path):
    with pytest.raises(ValueError,match='separate'):
        installer.build(tmp_path/'build/version',tmp_path/'build')
    with pytest.raises(ValueError,match='seed'):
        installer.build(tmp_path/'version',tmp_path/'seed/build',seed=tmp_path/'seed')


def test_relate_uses_dynamic_scheme_and_preserves_configured_cgmlst(tmp_path, monkeypatch):
    installed(tmp_path/'mlst_db')
    assemblies=tmp_path/'assemblies'; assemblies.mkdir()
    (assemblies/'S1.fasta').write_text('>contig\nACGT\n')
    sheet=tmp_path/'sheet.tsv';sheet.write_text('S1\tS1\tbrain\n')
    output=tmp_path/'result'
    monkeypatch.setattr(typing,'_check_tool',lambda *a,**kw:True)
    calls=[]
    def mlst(*args,**kwargs):
        calls.append(kwargs['scheme']); return 'ST1'
    monkeypatch.setattr(typing,'run_mlst',mlst)
    args=['--name','test','--run-name','test','--accession','test','--organism','Bacillus subtilis',
          '--assemblies-dir',str(assemblies),'--sample-sheet',str(sheet),'--outdir',str(output),
          '--db-dir',str(tmp_path),'--no-ani','--no-ska','--no-pdf']
    assert typing.main(args)==0
    assert calls==['bsubtilis_test']
    result=json.loads((output/'analysis_result.json').read_text())
    assert result['parameters']['mlst_selection']['scheme']=='bsubtilis_test'
    assert result['parameters']['cgmlst'] is False
    calls.clear()
    assert typing.main(args+['--resume'])==0
    assert not calls
    # Even if both schemes return ST1, an ST belongs to its chosen scheme.
    installed(tmp_path/'mlst_db',scheme='bsubtilis_alternative')
    assert typing.main(args+['--resume'])==0
    assert calls==['bsubtilis_alternative']
