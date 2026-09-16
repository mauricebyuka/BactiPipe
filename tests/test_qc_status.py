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
        assert interpret(expected, hit, ['Pass']*4)['overall'] == 'Review'
        assert interpret(expected, hit, ['Pass','Fail','Pass','Pass'])['overall'] == 'Fail'


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
