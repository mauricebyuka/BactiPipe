import csv
from pathlib import Path

from bactipipe.genome_sizes import GenomeSizeCatalog, coverage_and_verdict
from bactipipe.scripts.qc_plan import (
    bounded_worker_count,
    flye_size_arguments,
    illumina_minimum_coverage,
    nanopore_qc_failure_reason,
    overall_quality_verdict,
    write_execution_plan,
)


def test_controlled_names_resolve_case_and_whitespace_without_guessing():
    catalog = GenomeSizeCatalog.from_package()

    normal = catalog.resolve("Escherichia coli")
    harmless_variant = catalog.resolve("  ESCHERICHIA   COLI ")
    unsupported = catalog.resolve("Salmonella imaginary")

    assert normal.status == "known"
    assert harmless_variant == normal.__class__(
        entered_name="ESCHERICHIA   COLI",
        status="known",
        canonical_name=normal.canonical_name,
        taxid=normal.taxid,
        genome_size_bp=normal.genome_size_bp,
        reference_accession=normal.reference_accession,
        source=normal.source,
    )
    assert unsupported.status == "unsupported"
    assert unsupported.genome_size_bp is None


def test_unknown_and_lambda_have_explicit_distinct_behavior():
    catalog = GenomeSizeCatalog.from_package()

    unknown = catalog.resolve("unknown")
    lambda_control = catalog.resolve("lambda")

    assert unknown.status == "unknown"
    assert unknown.genome_size_bp is None
    assert flye_size_arguments(unknown.genome_size_bp) == []
    assert lambda_control.genome_size_bp == 48_502
    assert flye_size_arguments(lambda_control.genome_size_bp) == [
        "--asm-coverage",
        "50",
        "--g",
        "48502",
    ]


def test_coverage_unavailable_is_not_a_pass():
    catalog = GenomeSizeCatalog.from_package()
    assert coverage_and_verdict(10_000_000, catalog.resolve("Unknown"), 50) == (
        None,
        "N/A",
    )
    coverage, verdict = coverage_and_verdict(
        100_000_000, catalog.resolve("Escherichia coli"), 50
    )
    assert coverage is not None
    assert verdict == "Fail"


def test_worker_count_never_oversubscribes_explicit_cpu_allocation():
    assert bounded_worker_count(
        sample_count=12,
        total_cpus=24,
        cpus_per_sample=8,
        memory_worker_limit=10,
    ) == 3
    assert bounded_worker_count(
        sample_count=2,
        total_cpus=24,
        cpus_per_sample=8,
        memory_worker_limit=10,
    ) == 2


def test_nanopore_qc_gate_does_not_reject_low_estimated_depth():
    metrics = "Trimmed Total Bases: 100,000,000\nTrimmed Mean Read Quality: 16.2\n"
    assert nanopore_qc_failure_reason(
        metrics,
        genome_size_bp=1_000_000,
        minimum_quality=15,
        minimum_coverage=50,
    ) is None
    assert nanopore_qc_failure_reason(
        metrics,
        genome_size_bp=5_000_000,
        minimum_quality=15,
        minimum_coverage=50,
    ) is None
    assert nanopore_qc_failure_reason(
        metrics,
        genome_size_bp=None,
        minimum_quality=15,
        minimum_coverage=50,
    ) is None


def test_nanopore_qc_still_rejects_empty_low_quality_and_incomplete_reads():
    for metrics, reason in [
        ("Trimmed Total Bases: 0\nTrimmed Mean Read Quality: 16.2\n", "no sequence bases"),
        ("Trimmed Total Bases: 100\nTrimmed Mean Read Quality: 9.9\n", "mean read quality"),
        ("Trimmed Total Bases: 100\n", "complete QC outputs"),
    ]:
        assert reason in nanopore_qc_failure_reason(
            metrics, genome_size_bp=None, minimum_quality=15, minimum_coverage=50
        )


def test_illumina_coverage_threshold_uses_identified_species():
    assert illumina_minimum_coverage(30, "Salmonella enterica") == 30
    assert illumina_minimum_coverage(30, "Escherichia coli") == 40
    assert illumina_minimum_coverage(28, "Escherichia coli") == 60


def test_unresolved_coverage_cannot_produce_an_overall_pass():
    assert overall_quality_verdict(
        quality_verdict="Pass",
        coverage_verdict="N/A",
        taxonomy_verdict="N/A",
        completeness_verdict="Pass",
        contamination_verdict="Pass",
    ) == "Review"
    assert overall_quality_verdict(
        quality_verdict="Pass",
        coverage_verdict="Pass",
        taxonomy_verdict="N/A",
        completeness_verdict="Pass",
        contamination_verdict="Pass",
    ) == "Pass"


def test_execution_plan_records_per_sample_sizes_and_unknown_mode(tmp_path: Path):
    catalog = GenomeSizeCatalog.from_package()
    path = tmp_path / "sample_execution_plan.tsv"
    write_execution_plan(
        path,
        [
            ("lambda", "Lambda", catalog.resolve("Lambda")),
            ("ecoli", "Escherichia coli", catalog.resolve("Escherichia coli")),
            ("mystery", "Unknown", catalog.resolve("Unknown")),
        ],
        platform="Nanopore",
        min_quality=15,
        min_read_length=1000,
        min_coverage=50,
        desired_coverage=300,
        catalog_checksum=catalog.checksum,
    )
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, dialect="excel-tab"))
    assert [row["sample"] for row in rows] == ["lambda", "ecoli", "mystery"]
    assert rows[0]["genome_size_bp"] == "48502"
    assert rows[0]["qc_target_bases"] == str(48_502 * 300)
    assert rows[1]["genome_size_bp"] != rows[0]["genome_size_bp"]
    assert rows[2]["genome_size_bp"] == "N/A"
    assert rows[2]["assembly_size_mode"] == "size_independent"
    assert rows[2]["catalog_sha256"] == catalog.checksum


def test_illumina_execution_plan_is_explicitly_size_independent(tmp_path: Path):
    catalog = GenomeSizeCatalog.from_package()
    path = tmp_path / "sample_execution_plan.tsv"
    write_execution_plan(
        path,
        [("ecoli", "Escherichia coli", catalog.resolve("Escherichia coli"))],
        platform="Illumina",
        min_quality=28,
        min_read_length=None,
        min_coverage="quality_dependent_30_40_50_plus_ecoli_10",
        catalog_checksum=catalog.checksum,
        assembly_uses_genome_size=False,
    )
    with path.open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, dialect="excel-tab"))
    assert row["genome_size_bp"] != "N/A"
    assert row["assembly_size_mode"] == "size_independent"
    assert row["minimum_coverage_depth"] == (
        "quality_dependent_30_40_50_plus_ecoli_10"
    )
