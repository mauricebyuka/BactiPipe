import pytest
from bactipipe.genome_sizes import GenomeSizeCatalog
from bactipipe.organism_names import expected_organism_name
from bactipipe.qc_status import interpret
from bactipipe.species_identification import _taxonomic_confirmation


@pytest.mark.parametrize('entered,canonical', [
    ('Salmonella enterica Typhimurium', 'Salmonella enterica'),
    ('Salmonella enterica Dublin', 'Salmonella enterica'),
    ('Salmonella enterica Typhimurium strain LT2', 'Salmonella enterica'),
    ('Salmonella enterica serovar: Reading', 'Salmonella enterica'),
    ('Salmonella typhimurium', 'Salmonella enterica'),
    (' ESCHERICHIA  COLI strain K-12 ', 'Escherichia coli'),
    ('Escherichia coli K-12', 'Escherichia coli'),
    ('Escherichia coli O157:H7', 'Escherichia coli'),
    ('Staphylococcus aureus ATCC 25923', 'Staphylococcus aureus'),
    ('Staphylococcus aureus (isolate LabA)', 'Staphylococcus aureus'),
    ('Streptococcus equi ssp. equi strain 4047', 'Streptococcus equi subsp. equi'),
    ('Streptococcus equi subsp zooepidemicus isolate LabA', 'Streptococcus equi subsp. zooepidemicus'),
])
def test_shared_resolution_preserves_entry_and_uses_canonical_size(entered, canonical):
    catalog = GenomeSizeCatalog.from_package()
    result = catalog.resolve(entered)
    base = catalog.resolve(canonical)
    assert result.is_known
    assert result.entered_name == entered.strip()
    assert result.canonical_name == canonical
    assert result.genome_size_bp == base.genome_size_bp
    assert result.taxid == base.taxid
    status = interpret(entered, canonical, ['Pass'] * 4)
    assert status['overall'] == 'Pass'
    assert status['identification'] == 'matched'
    assert not status['input']['warning']
    assert _taxonomic_confirmation(entered, canonical) == 'Pass'


@pytest.mark.parametrize('entered', [
    'Staphylococcus epidermis', 'Campylobacter', 'Escherichia coli mysterious',
    'Escherichia coli / Shigella', 'Escherichia coli complex',
    'Escherichia coli or Klebsiella pneumoniae', 'Escherichia coli 2',
    'Streptococcus equi', 'Streptococcus equi subsp. invented strain 123',
    'Unknown strain K12',
])
def test_no_fuzzy_species_guessing_or_subspecies_fallback(entered):
    result = GenomeSizeCatalog.from_package().resolve(entered)
    assert result.is_unknown
    assert result.genome_size_bp is None


@pytest.mark.parametrize('expected,observed', [
    ('Salmonella enterica Typhimurium', 'Salmonella bongori'),
    ('Streptococcus equi subsp. equi strain 4047', 'Streptococcus equi subsp. zooepidemicus'),
    ('Streptococcus equi subsp. equi strain 4047', 'Streptococcus equi'),
])
def test_distinct_taxa_remain_distinct(expected, observed):
    assert _taxonomic_confirmation(expected, observed) == 'Fail'
    assert interpret(expected, observed, ['Pass'] * 4)['overall'] == 'Review'


def test_identification_is_not_limited_to_qc_catalog():
    name = 'Novelgenus novelspecies strain LAB123'
    assert expected_organism_name(name) == 'novelgenus novelspecies'
    assert not GenomeSizeCatalog.from_package().resolve(name).is_known
    assert _taxonomic_confirmation(name, 'Novelgenus novelspecies') == 'Pass'
