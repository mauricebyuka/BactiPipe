from bactipipe.qc_status import interpret, organism_input
from bactipipe.genome_sizes import GenomeSizeCatalog


def test_missing_and_unrecognized_inputs_use_unknown_without_size():
    catalog = GenomeSizeCatalog.from_package()
    for entered in ('', 'Campylobacter', 'novel bacterium'):
        result = catalog.resolve(entered)
        assert result.is_unknown and result.genome_size_bp is None
        assert result.canonical_name == 'Unknown'
        assert organism_input(entered)['warning']
    assert not organism_input('Unknown')['warning']


def test_review_does_not_mask_real_qc_failure():
    for expected, hit in [('Unknown', 'Closest: Escherichia coli'), ('Campylobacter', 'Escherichia coli'),
                          ('Salmonella enterica', 'Escherichia coli')]:
        assert interpret(expected, hit, ['Pass']*4)['overall'] == ('Review' if hit.startswith('Closest:') else 'Pass')
        assert interpret(expected, hit, ['Pass','Fail','Pass','Pass'])['overall'] == 'Fail'


def test_unexpected_confident_identity_is_pass_with_annotation():
    result = interpret('Salmonella enterica', 'Escherichia coli', ['Pass']*4)
    assert result['overall'] == 'Pass'
    assert result['identification'] == 'mismatch'
    for expected in ('', 'Typographical error'):
        result = interpret(expected, 'Escherichia coli', ['Pass']*4)
        assert result['overall'] == 'Pass' and result['input']['warning']


def test_intentional_unknown_and_matched_pass():
    assert interpret('Unknown', 'Escherichia coli', ['Pass']*4)['overall'] == 'Pass'
    assert interpret('Escherichia coli', 'Escherichia coli', ['Pass']*4)['overall'] == 'Pass'


def test_unavailable_qc_and_unresolved_identity_are_review():
    r = interpret('Unknown', 'No confident identification', ['Pass']*4)
    assert r['overall'] == 'Review' and r['quality'] == 'Pass'
    r = interpret('Unknown', 'Escherichia coli', ['Pass','N/A','Pass','Pass'])
    assert r['overall'] == 'Review' and r['quality'] == 'Not assessed'
    assert interpret('Unknown', 'Analysis failed', ['Pass']*4)['overall'] == 'Fail'


def test_contamination_precedes_identity_review():
    assert interpret('Salmonella enterica', 'Escherichia coli', ['Pass']*4,
                     isolate_qc_failure=True)['overall'] == 'Fail'


def test_nonfinite_measurements_remain_unassessed():
    from bactipipe.qc_status import threshold_verdict
    for value in [None, 'N/A', float('nan'), float('inf')]:
        assert threshold_verdict(value, 90) == 'N/A'


def test_summary_mismatch_passes_but_contamination_still_fails():
    import csv
    import tempfile
    from pathlib import Path
    from unittest.mock import patch
    from bactipipe.scripts.process_data import make_summary
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        source = root / 'input.tsv'
        target = root / 'summary.tsv'
        source.write_text('Sample\tQ\tqc\tExpected\tHit\tMatch\tDepth\tMin\tcov\ttax\n'
            'S1\t30\tPass\tSalmonella enterica\tEscherichia coli\t99\t100X\t50\tPass\tFail\n'
            'S2\t30\tPass\tSalmonella enterica\tEscherichia coli\t99\t100X\t50\tPass\tContamination\n'
            'S3\t30\tPass\tTypographical error\tEscherichia coli\t99\t100X\t50\tPass\tFail\n'
            'S4\t30\tPass\tUnknown\tEscherichia coli\t99\t100X\t50\tPass\tMixedSignalReview\n')
        with patch('bactipipe.scripts.process_data.render_pdf_from_tsv'), patch('bactipipe.scripts.process_data.simple_print'):
            make_summary(str(target), str(source), [], {'S1':['99','0'], 'S2':['99','0'], 'S3':['99','0'], 'S4':['99','0']}, logfile=str(root/'log'))
        lines = target.read_text().splitlines()
        start = next(i for i, line in enumerate(lines) if line.startswith('Sample ID\t'))
        rows = list(csv.DictReader(lines[start:], delimiter='\t'))
        assert [row['Overall Quality'] for row in rows] == ['Pass', 'Fail', 'Pass', 'Review']
        assert rows[3]['Isolate QC'] == 'Review'
        assert rows[3]['Identified organism'] == 'Escherichia coli'
