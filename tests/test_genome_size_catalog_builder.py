import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "build_controlled_genome_size_catalog.py"
SPEC = importlib.util.spec_from_file_location("genome_size_catalog_builder", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_chromosome_catalog_excludes_plasmids_and_uses_species_median():
    controlled = [
        {
            "taxid": "562",
            "legacy_size_bp": 5_000_000,
            "controlled_accession": "NZ_TEST",
            "canonical_name": "Escherichia coli",
        }
    ]
    accession_taxids = {"GCF_A": {"562"}, "GCF_B": {"562"}}
    sequence_reports = [
        {
            "assembly_accession": "GCF_A",
            "assigned_molecule_location_type": "Chromosome",
            "length": 4_700_000,
        },
        {
            "assembly_accession": "GCF_A",
            "assigned_molecule_location_type": "Plasmid",
            "length": 200_000,
        },
        {
            "assembly_accession": "GCF_B",
            "assigned_molecule_location_type": "Chromosome",
            "length": 5_100_000,
        },
    ]
    row = MODULE.summarize_catalog(
        controlled, accession_taxids, sequence_reports
    )[0]
    assert row["expected_size_bp"] == 4_900_000
    assert row["assembly_count"] == 2
    assert row["method"] == "median_refseq_chromosome_length"


def test_missing_chromosome_report_is_an_explicit_fallback():
    controlled = [
        {
            "taxid": "1",
            "legacy_size_bp": 1_234,
            "controlled_accession": "NZ_TEST",
            "canonical_name": "Example bacterium",
        }
    ]
    row = MODULE.summarize_catalog(controlled, {}, [])[0]
    assert row["expected_size_bp"] == 1_234
    assert row["assembly_count"] == 0
    assert row["method"] == "controlled_reference_fallback_no_chromosome_report"
