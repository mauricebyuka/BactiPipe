from __future__ import annotations

import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from bactipipe.species_identification import (
    AniHit,
    ContaminationEvidence,
    IdentificationPolicy,
    ReferenceCatalog,
    ReferenceRecord,
    TaxonomicSignal,
    assess_mixed_signals,
    checkm_contamination_evidence,
    default_kraken_confidence,
    default_mixed_signal_policy,
    decide_identification,
    parse_kmerfinder_output,
    parse_skani_output,
    write_identification_outputs,
)
from bactipipe.species_workflow import legacy_summary_fields
from bactipipe.species_workflow import identify_assembly


def hit(
    taxid: str,
    species: str,
    ani: float,
    af_ref: float,
    af_query: float,
    accession: str | None = None,
    reporting_group: str = "",
) -> AniHit:
    return AniHit(
        accession=accession or f"GCF_{taxid}",
        species_taxid=taxid,
        species_name=species,
        reference_path=f"{accession or taxid}.fna",
        ani=ani,
        aligned_fraction_reference=af_ref,
        aligned_fraction_query=af_query,
        reporting_group=reporting_group,
    )


class SpeciesDecisionTests(unittest.TestCase):
    def test_gtdb_ani_and_either_af_boundaries_for_both_platforms(self):
        for platform in ('Illumina', 'Nanopore'):
            for ani, ref, query, expected in (
                (95, 50, 40, 'confirmed_species'),
                (95, 40, 50, 'confirmed_species'),
                (95.44, 71.42, 62.49, 'confirmed_species'),
                (94.999, 90, 90, 'no_confident_identification'),
                (99, 49.999, 49.999, 'no_confident_identification'),
            ):
                with self.subTest(platform=platform, ani=ani, ref=ref, query=query):
                    result = decide_identification(sample_id='S', platform=platform,
                        ani_hits=[hit('44008', 'Enterococcus cecorum', ani, ref, query)])
                    self.assertEqual(result.status, expected)
                    self.assertEqual(result.identification_thresholds['minimum_aligned_fraction'], 50)
                    self.assertEqual(result.identification_thresholds['aligned_fraction_mode'], 'either')

    def test_legacy_bidirectional_mode_and_invalid_mode(self):
        with self.assertRaises(ValueError):
            IdentificationPolicy(aligned_fraction_mode='average')
        result = decide_identification(sample_id='S',
            policy=IdentificationPolicy(minimum_aligned_fraction=65, aligned_fraction_mode='both'),
            ani_hits=[hit('1', 'Species alpha', 99.2, 6.7, 95)])
        self.assertEqual(result.status, 'no_confident_identification')

    def test_confident_species_reports_ani_as_match(self):
        result = decide_identification(
            sample_id="S1",
            expected_organism="Mannheimia haemolytica",
            ani_hits=[
                hit("75985", "Mannheimia haemolytica", 98.7, 94.0, 92.0),
                hit("85404", "Mannheimia varigena", 92.2, 88.0, 84.0),
            ],
        )
        self.assertEqual(result.status, "confirmed_species")
        self.assertEqual(result.reported_organism, "Mannheimia haemolytica")
        self.assertEqual(result.match_percent, 98.7)
        self.assertEqual(result.taxonomic_confirmation, "Pass")

    def test_salmonella_serovar_expectation_is_compared_at_species_level(self):
        result = decide_identification(
            sample_id="S1",
            expected_organism="Salmonella enterica Dublin",
            ani_hits=[hit("28901", "Salmonella enterica", 99.8, 96.0, 95.0)],
        )
        self.assertEqual(result.status, "confirmed_species")
        self.assertEqual(result.reported_organism, "Salmonella enterica")
        self.assertEqual(result.taxonomic_confirmation, "Pass")

    def test_historical_salmonella_typhimurium_label_is_species_concordant(self):
        result = decide_identification(
            sample_id="S1",
            expected_organism="Salmonella typhimurium",
            ani_hits=[hit("28901", "Salmonella enterica", 99.8, 96.0, 95.0)],
        )
        self.assertEqual(result.taxonomic_confirmation, "Pass")

    def test_species_normalization_does_not_merge_salmonella_bongori(self):
        result = decide_identification(
            sample_id="S1",
            expected_organism="Salmonella enterica Dublin",
            ani_hits=[hit("54736", "Salmonella bongori", 99.8, 96.0, 95.0)],
        )
        self.assertEqual(result.taxonomic_confirmation, "Fail")

    def test_streptococcus_equi_subspecies_requires_distinct_reported_name(self):
        result = decide_identification(
            sample_id="S1",
            expected_organism="Streptococcus equi subsp. equi",
            ani_hits=[hit("1336", "Streptococcus equi", 99.8, 96.0, 95.0)],
        )
        self.assertEqual(result.taxonomic_confirmation, "Fail")

        resolved = decide_identification(
            sample_id="S1",
            expected_organism="Streptococcus equi subsp. equi",
            ani_hits=[
                hit(
                    "1336",
                    "Streptococcus equi subsp. equi",
                    99.8,
                    96.0,
                    95.0,
                )
            ],
        )
        self.assertEqual(resolved.taxonomic_confirmation, "Pass")

    def test_aligned_fraction_advantage_resolves_close_ani_candidates(self):
        result = decide_identification(
            sample_id="S1",
            expected_organism="Streptococcus equi subsp. zooepidemicus",
            ani_hits=[
                hit(
                    "40041",
                    "Streptococcus equi subsp. zooepidemicus",
                    97.58,
                    88.0,
                    90.0,
                ),
                hit(
                    "1336",
                    "Streptococcus equi subsp. equi",
                    97.31,
                    80.0,
                    88.0,
                ),
            ],
        )
        self.assertEqual(result.status, "confirmed_species")
        self.assertEqual(
            result.reported_organism,
            "Streptococcus equi subsp. zooepidemicus",
        )
        self.assertEqual(result.taxonomic_confirmation, "Pass")

    def test_low_genome_fraction_is_not_a_positive_identification(self):
        result = decide_identification(
            sample_id="P2631472",
            expected_organism="Unknown",
            ani_hits=[hit("75985", "Mannheimia haemolytica", 99.2, 6.7, 7.3)],
        )
        self.assertEqual(result.status, "no_confident_identification")
        self.assertEqual(result.reported_organism, "No confident identification")
        self.assertIsNone(result.match_percent)
        self.assertEqual(result.closest_species, ())

    def test_close_passing_species_are_not_forced(self):
        result = decide_identification(
            sample_id="S1",
            ani_hits=[
                hit("1", "Species alpha", 97.10, 90.0, 90.0),
                hit("2", "Species beta", 96.75, 91.0, 89.0),
            ],
        )
        self.assertEqual(result.status, "ambiguous_species")
        self.assertIsNone(result.match_percent)

    def test_best_reference_is_selected_within_species(self):
        result = decide_identification(
            sample_id="S1",
            ani_hits=[
                hit("1", "Species alpha", 96.0, 90.0, 90.0, "GCF_A"),
                hit("1", "Species alpha", 99.0, 70.0, 70.0, "GCF_B"),
                hit("2", "Species beta", 94.0, 95.0, 95.0, "GCF_C"),
            ],
        )
        self.assertEqual(result.status, "confirmed_species")
        self.assertEqual(result.match_percent, 99.0)
        self.assertEqual(result.closest_species[0].accession, "GCF_B")

    def test_reporting_group_prevents_false_ani_ambiguity(self):
        group = "Escherichia coli/Shigella complex"
        result = decide_identification(
            sample_id="C04",
            expected_organism="Escherichia coli",
            ani_hits=[
                hit("562", "Escherichia coli", 98.70, 96.0, 96.0, reporting_group=group),
                hit("623", "Shigella flexneri", 99.75, 98.0, 98.0, reporting_group=group),
            ],
            reporting_group_members={
                group: ("Escherichia coli", "Shigella flexneri")
            },
        )
        self.assertEqual(result.status, "confirmed_species")
        self.assertEqual(result.reported_organism, "Escherichia coli")
        self.assertEqual(result.match_percent, 98.70)
        self.assertEqual(result.taxonomic_confirmation, "Pass")
        self.assertEqual(legacy_summary_fields(result).coverage_organism, "Escherichia coli")
        self.assertEqual(result.closest_species[0].reporting_group, group)

    def test_reporting_group_is_supporting_evidence_for_unknown_sample(self):
        group = "Escherichia coli/Shigella complex"
        result = decide_identification(
            sample_id="unknown",
            expected_organism="Unknown",
            ani_hits=[
                hit("562", "Escherichia coli", 98.70, 96.0, 96.0, reporting_group=group),
                hit("623", "Shigella flexneri", 99.75, 98.0, 98.0, reporting_group=group),
            ],
            reporting_group_members={
                group: ("Escherichia coli", "Shigella flexneri")
            },
        )
        self.assertEqual(result.reported_organism, "Shigella flexneri")
        self.assertNotEqual(result.reported_organism, group)
        self.assertEqual(result.closest_species[0].reporting_group, group)

    def test_reporting_group_collapses_mixed_signal_evidence(self):
        group = "Escherichia coli/Shigella complex"
        signals = [
            TaxonomicSignal("562", "Escherichia coli", 7000, 55.0, 90.0, 6000),
            TaxonomicSignal("623", "Shigella flexneri", 5000, 40.0, 85.0, 4500),
        ]
        without_group = assess_mixed_signals(signals)
        with_group = assess_mixed_signals(
            signals,
            reporting_groups={"562": group, "623": group},
        )
        self.assertEqual(without_group.status, "strong_secondary")
        self.assertEqual(with_group.status, "single_signal")

    def test_reporting_group_does_not_hide_unrelated_secondary_species(self):
        group = "Escherichia coli/Shigella complex"
        assessment = assess_mixed_signals(
            [
                TaxonomicSignal("562", "Escherichia coli", 7000, 50.0, 90.0, 6000),
                TaxonomicSignal("623", "Shigella flexneri", 5000, 35.0, 85.0, 4500),
                TaxonomicSignal("28901", "Salmonella enterica", 4000, 15.0, 80.0, 3500),
            ],
            reporting_groups={"562": group, "623": group},
        )
        self.assertEqual(assessment.status, "strong_secondary")
        self.assertEqual(
            assessment.qualifying_secondary[0].species_name,
            "Salmonella enterica",
        )

    def test_shared_ncbi_species_parent_groups_subspecies_only_for_mixed_signal(self):
        records = [
            ReferenceRecord(
                accession="GCA_EQUI",
                species_taxid="1336",
                species_name="Streptococcus equi subsp. equi",
                fasta_path="equi.fna",
                source="GenBank",
                reference_role="type",
                assembly_level="Complete Genome",
                sequence_origin="isolate",
                type_material=True,
                canonical_for_tree=True,
                callable_for_identification=True,
                ncbi_taxonomy_name="Streptococcus equi",
            ),
            ReferenceRecord(
                accession="GCF_ZOO",
                species_taxid="40041",
                species_name="Streptococcus equi subsp. zooepidemicus",
                fasta_path="zoo.fna",
                source="RefSeq",
                reference_role="type",
                assembly_level="Complete Genome",
                sequence_origin="isolate",
                type_material=True,
                canonical_for_tree=True,
                callable_for_identification=True,
                ncbi_taxonomy_name="Streptococcus equi",
            ),
        ]
        catalog = ReferenceCatalog(records, root=Path("."))
        self.assertEqual(
            catalog.mixed_signal_groups_by_taxid["1336"],
            catalog.mixed_signal_groups_by_taxid["40041"],
        )

        assessment = assess_mixed_signals(
            [
                TaxonomicSignal(
                    "1336", records[0].species_name, 638373, 40.3, 60.2, 435403
                ),
                TaxonomicSignal(
                    "40041", records[1].species_name, 895852, 56.6, 76.3, 112128
                ),
                TaxonomicSignal("28901", "Salmonella enterica", 20, 0.1, 0.2, 10),
            ],
            reporting_groups=catalog.mixed_signal_groups_by_taxid,
        )
        self.assertEqual(assessment.status, "single_signal")

        result = decide_identification(
            sample_id="S1",
            expected_organism=records[1].species_name,
            ani_hits=[
                hit("40041", records[1].species_name, 97.5, 87.0, 82.0),
                hit("1336", records[0].species_name, 97.3, 80.0, 78.0),
            ],
            mixed_signal=assessment,
        )
        self.assertEqual(result.status, "confirmed_species")
        self.assertEqual(result.reported_organism, records[1].species_name)

    def test_strong_secondary_requires_cohort_review_before_contamination_call(self):
        mixed = assess_mixed_signals(
            [
                TaxonomicSignal("1", "Species alpha", 10000, 60.0, 90.0, 9000),
                TaxonomicSignal("2", "Species beta", 4000, 25.0, 75.0, 3500),
            ]
        )
        self.assertEqual(mixed.status, "strong_secondary")
        preliminary = decide_identification(
            sample_id="S1",
            ani_hits=[hit("1", "Species alpha", 99.0, 95.0, 95.0)],
            mixed_signal=mixed,
        )
        self.assertEqual(preliminary.status, "mixed_signal_review")
        result = decide_identification(
            sample_id="S1",
            ani_hits=[hit("1", "Species alpha", 99.0, 95.0, 95.0)],
            mixed_signal=mixed,
            cohort_review_complete=True,
        )
        self.assertEqual(result.status, "suspected_culture_contamination")
        self.assertEqual(result.reported_organism, "Species alpha")
        self.assertEqual(result.match_percent, 99.0)

    def test_illumina_mapping_alone_cannot_declare_culture_contamination(self):
        mixed = assess_mixed_signals(
            [
                TaxonomicSignal("1", "Species alpha", 10000, 60.0, 90.0, 9000),
                TaxonomicSignal("2", "Species beta", 4000, 25.0, 75.0, 3500),
            ]
        )
        result = decide_identification(
            sample_id="S1",
            expected_organism="Species alpha",
            ani_hits=[hit("1", "Species alpha", 99.0, 95.0, 95.0)],
            mixed_signal=mixed,
            cohort_review_complete=True,
            platform="Illumina",
            contamination_evidence=checkm_contamination_evidence(1.0),
        )
        self.assertEqual(result.status, "mixed_signal_review")
        self.assertIn("mapping alone is insufficient", result.reason)

    def test_illumina_checkm_corroborates_culture_contamination(self):
        mixed = assess_mixed_signals(
            [
                TaxonomicSignal("1", "Species alpha", 10000, 60.0, 90.0, 9000),
                TaxonomicSignal("2", "Species beta", 4000, 25.0, 75.0, 3500),
            ]
        )
        result = decide_identification(
            sample_id="S1",
            expected_organism="Species alpha",
            ani_hits=[hit("1", "Species alpha", 99.0, 95.0, 95.0)],
            mixed_signal=mixed,
            cohort_review_complete=True,
            platform="Illumina",
            contamination_evidence=checkm_contamination_evidence(99.84),
        )
        self.assertEqual(result.status, "suspected_culture_contamination")
        self.assertEqual(result.reported_organism, "Species alpha")
        self.assertEqual(result.taxonomic_confirmation, "Pass")
        self.assertEqual(
            result.contamination_evidence.status, "supports_contamination"
        )

    def test_missing_checkm_keeps_illumina_secondary_signal_in_review(self):
        evidence = checkm_contamination_evidence("N/A")
        self.assertEqual(evidence, ContaminationEvidence())

    def test_broad_but_incomplete_secondary_remains_review_only(self):
        assessment = assess_mixed_signals(
            [
                TaxonomicSignal("1", "Species alpha", 10000, 55.0, 51.0, 9000),
                TaxonomicSignal("2", "Species beta", 8000, 40.0, 37.0, 7000),
            ],
            policy=default_mixed_signal_policy("Illumina"),
        )
        self.assertEqual(assessment.status, "minor_secondary")
        result = decide_identification(
            sample_id="S1",
            ani_hits=[hit("1", "Species alpha", 90.0, 45.0, 40.0)],
            mixed_signal=assessment,
            cohort_review_complete=True,
        )
        self.assertEqual(result.status, "mixed_signal_review")

    def test_shared_hit_without_breadth_cannot_declare_mixture(self):
        assessment = assess_mixed_signals(
            [
                TaxonomicSignal("1", "Species alpha", 10000, 90.0, 95.0, 9000),
                TaxonomicSignal("2", "Species beta", 500, 5.0, None, 400),
            ]
        )
        self.assertEqual(assessment.status, "indeterminate")

    def test_tiny_unmapped_classifier_tail_does_not_force_review(self):
        assessment = assess_mixed_signals(
            [
                TaxonomicSignal("1", "Species alpha", 10000, 95.0, 95.0, 9000),
                TaxonomicSignal("2", "Species beta", 12, 0.1, None, 10),
            ]
        )
        self.assertEqual(assessment.status, "single_signal")

    def test_supported_minor_secondary_prevents_forced_species_call(self):
        assessment = assess_mixed_signals(
            [
                TaxonomicSignal("1", "Species alpha", 10000, 90.0, 95.0, 9000),
                TaxonomicSignal("2", "Species beta", 500, 5.0, 50.0, 400),
            ]
        )
        self.assertEqual(assessment.status, "minor_secondary")
        result = decide_identification(
            sample_id="S1",
            ani_hits=[hit("1", "Species alpha", 99.0, 95.0, 95.0)],
            mixed_signal=assessment,
        )
        self.assertEqual(result.status, "mixed_signal_review")
        self.assertIsNone(result.match_percent)

    def test_competitive_mapping_counts_override_classifier_counts(self):
        assessment = assess_mixed_signals(
            [
                TaxonomicSignal(
                    "1", "Species alpha", 10000, 70.0, 95.0, 9000, 9000, 99.0
                ),
                TaxonomicSignal(
                    "2", "Species beta", 4000, 25.0, 40.0, 3500, 20, 0.2
                ),
            ]
        )
        self.assertEqual(assessment.status, "single_signal")

    def test_unmapped_classifier_candidate_cannot_outrank_mapped_evidence(self):
        assessment = assess_mixed_signals(
            [
                TaxonomicSignal("9", "Spurious candidate", 50000, 50.0, None, 8000),
                TaxonomicSignal(
                    "1", "Species alpha", 20000, 40.0, 95.0, 9000, 18000, 100.0
                ),
            ]
        )
        self.assertEqual(assessment.status, "indeterminate")
        self.assertEqual(assessment.primary.species_name, "Species alpha")

    def test_illumina_mixed_signal_defaults_are_more_conservative(self):
        illumina = default_mixed_signal_policy("Illumina")
        nanopore = default_mixed_signal_policy("Nanopore")
        self.assertGreater(
            illumina.minimum_secondary_distinct_minimizers,
            nanopore.minimum_secondary_distinct_minimizers,
        )
        self.assertGreater(
            illumina.minimum_secondary_reference_breadth,
            nanopore.minimum_secondary_reference_breadth,
        )
        self.assertGreater(
            default_kraken_confidence("Illumina"),
            default_kraken_confidence("Nanopore"),
        )

    def test_legacy_summary_uses_ani_only_for_confirmed_call(self):
        result = decide_identification(
            sample_id="S1",
            ani_hits=[hit("1", "Species alpha", 99.0, 95.0, 93.0)],
        )
        fields = legacy_summary_fields(result)
        self.assertEqual(fields.identified_organism, "Species alpha")
        self.assertEqual(fields.match_percent, "99.00")
        self.assertEqual(fields.coverage_organism, "Species alpha")

    def test_legacy_summary_does_not_label_closest_ani_as_match(self):
        result = decide_identification(
            sample_id="S1",
            expected_organism="Species expected",
            ani_hits=[hit("1", "Species alpha", 93.0, 80.0, 75.0)],
        )
        fields = legacy_summary_fields(result)
        self.assertEqual(fields.identified_organism, "No confident identification")
        self.assertEqual(fields.match_percent, "N/A")
        self.assertEqual(fields.taxonomic_confirmation, "Fail")

    def test_unresolved_summary_does_not_replace_identity_with_supporting_evidence(self):
        from dataclasses import replace, asdict
        base = decide_identification(
            sample_id="S1", ani_hits=[hit("1", "Species alpha", 93.0, 80.0, 75.0)],
        )
        for status in ("mixed_signal_review", "mixed", "ambiguous_species",
                       "suspected_culture_contamination", "no_confident_identification"):
            result = replace(base, status=status)
            original = asdict(result)
            fields = legacy_summary_fields(result)
            self.assertEqual(fields.identified_organism, "No confident identification")
            self.assertEqual(fields.match_percent, "N/A")
            self.assertEqual(asdict(result), original)


class EvidenceParserTests(unittest.TestCase):
    def test_kmerfinder_coverage_is_preserved_as_coverage_not_identity(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            path = Path(temporary_dir) / "results.txt"
            path.write_text(
                "# Assembly\tNum\tScore\tExpected\tTemplate_length\tQuery_Coverage\t"
                "Template_Coverage\tDepth\ttot_query_Coverage\t"
                "tot_template_Coverage\ttot_depth\tq_value\tp_value\tAccession Number\t"
                "Description\tTAXID\tTaxonomy\tTAXID Species\tSpecies\n"
                "GCF_1\t1\t200\t0\t1000\t6.70\t95.00\t1\t6.70\t95.00\t1\t1\t1e-9\t"
                "NC_1\tdescription\t1\tBacteria\t75985\tMannheimia haemolytica\n"
            )
            evidence = parse_kmerfinder_output(path)
        self.assertEqual(len(evidence), 1)
        self.assertEqual(evidence[0].query_coverage, 6.7)
        self.assertEqual(evidence[0].template_coverage, 95.0)

    def test_manifest_controls_skani_taxonomy(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            fasta = root / "refs" / "GCF_1.fna"
            fasta.parent.mkdir()
            fasta.write_text(">ref\nACGT\n")
            checksum = hashlib.sha256(fasta.read_bytes()).hexdigest()
            manifest = root / "reference_manifest.tsv"
            manifest.write_text(
                "accession\tspecies_taxid\tspecies_name\tfasta_path\tsource\t"
                "reference_role\tassembly_level\tsequence_origin\ttype_material\t"
                "canonical_for_tree\tcallable_for_identification\t"
                "reporting_group\tsha256\n"
                f"GCF_1\t1\tSpecies alpha\trefs/GCF_1.fna\tRefSeq\ttype\tcomplete\t"
                f"isolate\tyes\tyes\tyes\t\t{checksum}\n"
            )
            catalog = ReferenceCatalog.from_tsv(manifest)
            catalog.validate_files()
            output = root / "skani.tsv"
            output.write_text(
                "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\t"
                "Align_fraction_query\tRef_name\tQuery_name\n"
                f"{fasta}\tquery.fna\t98.1\t90.0\t88.0\tref\tquery\n"
            )
            hits = parse_skani_output(output, catalog)
        self.assertEqual(hits[0].species_name, "Species alpha")
        self.assertEqual(hits[0].accession, "GCF_1")

    def test_manifest_resolves_reference_after_atomic_database_move(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            fasta = root / "GCF_000005845.2.fna"
            fasta.write_text(">ref\nACGT\n")
            checksum = hashlib.sha256(fasta.read_bytes()).hexdigest()
            manifest = root / "reference_manifest.tsv"
            manifest.write_text(
                "accession\tspecies_taxid\tspecies_name\tfasta_path\tsource\t"
                "reference_role\tassembly_level\tsequence_origin\ttype_material\t"
                "canonical_for_tree\tcallable_for_identification\treporting_group\tsha256\n"
                f"GCF_000005845.2\t562\tEscherichia coli\t{fasta}\tRefSeq\ttype\t"
                f"complete\tisolate\tyes\tyes\tyes\t\t{checksum}\n"
            )
            catalog = ReferenceCatalog.from_tsv(manifest)
            resolved = catalog.resolve(
                "/old/.candidate.building/references/GCF_000005845.2.fna"
            )
        self.assertEqual(resolved.species_taxid, "562")

    def test_writes_machine_readable_and_tabular_outputs(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            result = decide_identification(
                sample_id="S1",
                ani_hits=[hit("1", "Species alpha", 99.0, 90.0, 90.0)],
            )
            write_identification_outputs(result, temporary_dir)
            payload = json.loads(
                (Path(temporary_dir) / "species_identification.json").read_text()
            )
            self.assertEqual(payload["status"], "confirmed_species")
            self.assertEqual(payload["identification_thresholds"]["minimum_ani"], 95.0)
            with (Path(temporary_dir) / "closest_species.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["species_name"], "Species alpha")


class AtomicWorkflowPublicationTests(unittest.TestCase):
    def test_complete_result_replaces_stale_sample_evidence(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            final = Path(temporary_dir) / "S1"
            final.mkdir()
            (final / "stale.txt").write_text("old")
            expected = decide_identification(
                sample_id="S1",
                ani_hits=[hit("1", "Species alpha", 99.0, 90.0, 90.0)],
            )

            def fake_once(*, output_directory, **_kwargs):
                output = Path(output_directory)
                output.mkdir()
                (output / "fresh.txt").write_text("new")
                return expected

            with patch(
                "bactipipe.species_workflow._identify_assembly_once",
                side_effect=fake_once,
            ):
                observed = identify_assembly(output_directory=final)

            self.assertEqual(observed, expected)
            self.assertFalse((final / "stale.txt").exists())
            self.assertEqual((final / "fresh.txt").read_text(), "new")

    def test_failed_refresh_preserves_previous_complete_evidence(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            final = Path(temporary_dir) / "S1"
            final.mkdir()
            (final / "result.txt").write_text("previous")

            def fail_once(*, output_directory, **_kwargs):
                output = Path(output_directory)
                output.mkdir()
                (output / "partial.txt").write_text("partial")
                raise RuntimeError("test failure")

            with patch(
                "bactipipe.species_workflow._identify_assembly_once",
                side_effect=fail_once,
            ):
                with self.assertRaisesRegex(RuntimeError, "test failure"):
                    identify_assembly(output_directory=final)

            self.assertEqual((final / "result.txt").read_text(), "previous")
            self.assertFalse(any(final.parent.glob(".S1.building-*")))


if __name__ == "__main__":
    unittest.main()
