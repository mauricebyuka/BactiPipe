import json
import logging
from pathlib import Path

from bactipipe.scripts import type_genomes as typing


def test_cgmlst_timeout_is_unavailable_not_stale_success(tmp_path, monkeypatch):
    (tmp_path / 'SALMONELLA_summary.txt').write_text('cgST\n123\n')
    seen = {}
    def fake(cmd, **kwargs):
        seen.update(kwargs)
        return 124, 'partial evidence\nTimed out'
    monkeypatch.setattr(typing, '_run', fake)
    value = typing.run_cgmlst(tmp_path / 'S1.fasta', tmp_path, tmp_path, tmp_path,
                              logging.getLogger(__name__), 'SALMONELLA', timeout=17)
    assert value == 'NA' and seen['timeout'] == 17
    assert json.loads((tmp_path / 'cgmlst_status.json').read_text())['status'] == 'timed_out'
    assert 'partial evidence' in (tmp_path / 'cgmlst.stdout.txt').read_text()


def test_successful_cgmlst_is_unchanged(tmp_path, monkeypatch):
    (tmp_path / 'SALMONELLA_summary.txt').write_text('cgST\n123\n')
    monkeypatch.setattr(typing, '_run', lambda *a, **kw: (0, 'success'))
    assert typing.run_cgmlst(tmp_path / 'S1.fasta', tmp_path, tmp_path, tmp_path,
                             logging.getLogger(__name__), 'SALMONELLA') == '123'


def test_invalid_summary_records_failure(tmp_path, monkeypatch):
    (tmp_path / 'SALMONELLA_summary.txt').write_text('')
    monkeypatch.setattr(typing, '_run', lambda *a, **kw: (0, 'success'))
    assert typing.run_cgmlst(tmp_path / 'S1.fasta', tmp_path, tmp_path, tmp_path,
                             logging.getLogger(__name__), 'SALMONELLA') == 'NA'
    assert json.loads((tmp_path / 'cgmlst_status.json').read_text())['status'] == 'failed'


def test_timeout_cli_default_and_override():
    parser = typing._build_argparser()
    assert parser.parse_args([]).cgmlst_timeout == 900
    assert parser.parse_args(['--cgmlst-timeout', '1800']).cgmlst_timeout == 1800


def test_cohort_continues_after_timeout_and_resume_retries_only_unavailable(tmp_path, monkeypatch):
    for database in ('mlst_db', 'cgmlstfinder_db'):
        (tmp_path / database).mkdir()
    assemblies = tmp_path / 'assemblies'
    assemblies.mkdir()
    for sample in ('S1', 'S2'):
        (assemblies / f'{sample}.fasta').write_text('>contig\nACGTACGT\n')
    sheet = tmp_path / 'samples.tsv'
    sheet.write_text('S1\tFirst\tbrain\nS2\tSecond\tbrain\n')
    out = tmp_path / 'output'
    monkeypatch.setattr(typing, '_check_tool', lambda *a, **kw: True)
    monkeypatch.setattr(typing, 'serotype_dispatch', lambda *a, **kw: ('test', 'formula', 'Salmonella'))
    monkeypatch.setattr(typing, 'run_mlst', lambda *a, **kw: '1')
    called = []
    def tool(cmd, **kwargs):
        sample = Path(cmd[cmd.index('-i') + 1]).stem
        called.append(sample)
        if sample == 'S1':
            return 124, 'Timed out'
        folder = Path(cmd[cmd.index('-o') + 1])
        scheme = cmd[cmd.index('-s') + 1]
        (folder / f'{scheme}_summary.txt').write_text('cgST\n123\n')
        return 0, 'done'
    monkeypatch.setattr(typing, '_run', tool)
    args = ['--name', 'Tester', '--run-name', 'test', '--accession', 'cohort',
            '--organism', 'salmonella', '--assemblies-dir', str(assemblies),
            '--sample-sheet', str(sheet), '--outdir', str(out), '--db-dir', str(tmp_path),
            '--no-ani', '--no-ska', '--no-pdf']
    assert typing.main(args) == 0
    assert called == ['S1', 'S2']
    assert (out / 'cohort_relate.tsv').is_file()
    result = json.loads((out / 'analysis_result.json').read_text())
    assert result['status'] == 'completed'
    assert result['warnings'][0]['sample_id'] == 'S1'
    assert result['warnings'][0]['status'] == 'timed_out'
    assert 'S1\tNA' in (out / 'cohort_cgmlst.tsv').read_text()
    assert 'S2\t123' in (out / 'cohort_cgmlst.tsv').read_text()
    called.clear()
    assert typing.main(args + ['--resume']) == 0
    assert called == ['S1']
