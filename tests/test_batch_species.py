from __future__ import annotations

import csv
import hashlib
import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from bactipipe.batch_species import (
    _publish_reconciled_placement,
    classify_cross_sample_carryover,
    reconcile_batch_identifications,
    update_final_summary_verdicts,
    update_summary_identifications,
)
from bactipipe.phylogenetic_placement import PlacementResult
from bactipipe.kraken_evidence import write_taxonomic_signals
from bactipipe.species_identification import (
    TaxonomicSignal,
    decide_identification,
    write_identification_outputs,
)


def signal(taxid: str, name: str, reads: int, percent: float) -> TaxonomicSignal:
    return TaxonomicSignal(
        taxid,
        name,
        reads,
        percent,
        reference_breadth=80.0,
        distinct_minimizers=max(reads, 100),
        competitive_mapped_reads=reads,
        competitive_mapped_percent=percent,
    )


class CrossSampleClassificationTests(unittest.TestCase):
    def test_mapping_equivalent_subspecies_are_one_cross_sample_signal(self):
        group = "ncbi_species:streptococcus equi"
        evidence = classify_cross_sample_carryover(
            {
                "sample": [
                    signal("1336", "Streptococcus equi subsp. equi", 640_000, 40.0),
                    signal(
                        "40041",
                        "Streptococcus equi subsp. zooepidemicus",
                        900_000,
                        56.0,
                    ),
                ],
                "other": [
                    signal(
                        "40041",
                        "Streptococcus equi subsp. zooepidemicus",
                        1_000_000,
                        99.0,
                    )
                ],
            },
            signal_groups={"1336": group, "40041": group},
        )

        self.assertEqual(len(evidence["sample"]), 1)
        self.assertEqual(evidence["sample"][0].status, "source")
        self.assertEqual(evidence["sample"][0].reporting_group, "")

    def test_weak_secondary_with_strong_donor_is_possible_cross_talk(self):
        evidence = classify_cross_sample_carryover(
            {
                "salmonella_donor": [signal("28901", "Salmonella enterica", 100_000, 100.0)],
                "ecoli_recipient": [
                    signal("562", "Escherichia coli", 10_000, 95.0),
                    signal("28901", "Salmonella enterica", 500, 5.0),
                ],
            }
        )
        salmonella = next(
            item
            for item in evidence["ecoli_recipient"]
            if item.species_taxid == "28901"
        )
        self.assertEqual(salmonella.status, "possible_cross_talk")
        self.assertEqual(salmonella.source_sample, "salmonella_donor")
        self.assertAlmostEqual(salmonella.fraction_of_source, 0.005)

    def test_substantial_secondary_is_not_hidden_by_donor(self):
        evidence = classify_cross_sample_carryover(
            {
                "salmonella_donor": [signal("28901", "Salmonella enterica", 100_000, 100.0)],
                "mixed": [
                    signal("562", "Escherichia coli", 10_000, 70.0),
                    signal("28901", "Salmonella enterica", 4_000, 30.0),
                ],
            }
        )
        salmonella = next(
            item for item in evidence["mixed"] if item.species_taxid == "28901"
        )
        self.assertEqual(salmonella.status, "observed")

    def test_signal_without_cohort_donor_is_not_cross_talk(self):
        evidence = classify_cross_sample_carryover(
            {
                "sample": [
                    signal("562", "Escherichia coli", 10_000, 95.0),
                    signal("28901", "Salmonella enterica", 500, 5.0),
                ]
            }
        )
        salmonella = next(
            item for item in evidence["sample"] if item.species_taxid == "28901"
        )
        self.assertEqual(salmonella.status, "observed")

    def test_primary_donor_is_used_even_if_another_sample_has_more_secondary_reads(self):
        evidence = classify_cross_sample_carryover(
            {
                "salmonella_donor": [
                    signal("28901", "Salmonella enterica", 1_000, 100.0)
                ],
                "substantial_secondary": [
                    signal("562", "Escherichia coli", 20_000, 80.0),
                    signal("28901", "Salmonella enterica", 5_000, 20.0),
                ],
                "weak_recipient": [
                    signal("562", "Escherichia coli", 10_000, 99.5),
                    signal("28901", "Salmonella enterica", 50, 0.5),
                ],
            }
        )
        salmonella = next(
            item
            for item in evidence["weak_recipient"]
            if item.species_taxid == "28901"
        )
        self.assertEqual(salmonella.status, "possible_cross_talk")
        self.assertEqual(salmonella.source_sample, "salmonella_donor")


class BatchReconciliationTests(unittest.TestCase):
    def test_control_only_or_no_bacterial_results_writes_empty_audit(self):
        for expected in ({"Lambda": "Lambda"}, {"Lambda": "Lambda", "S1": "Salmonella enterica"}):
            with self.subTest(expected=expected), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary) / "species_identification"
                result = reconcile_batch_identifications(
                    species_output_root=root,
                    database_root=Path(temporary) / "database",
                    expected_organisms=expected,
                    platform="Nanopore",
                )
                self.assertEqual(result, {})
                audit = json.loads((root / "cross_sample_signal_review.json").read_text())
                self.assertEqual(audit["sample_count"], 0)
                self.assertEqual(len((root / "cross_sample_signal_review.tsv").read_text().splitlines()), 1)

    def test_unresolved_tree_is_published_after_cohort_decision(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample_dir = root / "species_identification" / "S1"
            sample_dir.mkdir(parents=True)
            assembly = root / "S1.fasta"
            assembly.write_text(">S1\nACGT\n")
            metadata = root / "metadata.tsv"
            metadata.write_text("accession\tgtdb_taxonomy\n")
            result = decide_identification(
                sample_id="S1",
                expected_organism="Unknown",
                ani_hits=(),
            )
            placement = PlacementResult(
                status="available",
                reason="placed",
                tree_newick="tree.nwk",
                tree_svg="tree.svg",
                closest_species=("Species alpha",),
            )
            with patch(
                "bactipipe.batch_species.run_bac120_placement",
                return_value=placement,
            ) as runner:
                payload = _publish_reconciled_placement(
                    sample_id="S1",
                    result=result,
                    sample_directory=sample_dir,
                    assembly=assembly,
                    gtdb_metadata_tsv=metadata,
                    gtdbtk_command_prefix=("gtdbtk",),
                    placement_threads=4,
                )
            self.assertEqual(payload["status"], "available")
            runner.assert_called_once()
            self.assertEqual(runner.call_args.kwargs["cpus"], 4)

    def test_mixed_tree_is_suppressed_after_cohort_decision(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample_dir = root / "species_identification" / "S1"
            stale = sample_dir / "phylogenetic_placement"
            stale.mkdir(parents=True)
            (stale / "closest_species_tree.svg").write_text("stale")
            result = decide_identification(
                sample_id="S1",
                expected_organism="Unknown",
                ani_hits=(),
            )
            result = replace(result, status="mixed_signal_review")
            payload = _publish_reconciled_placement(
                sample_id="S1",
                result=result,
                sample_directory=sample_dir,
                assembly=None,
                gtdb_metadata_tsv=root / "metadata.tsv",
                gtdbtk_command_prefix=("gtdbtk",),
                placement_threads=1,
            )
            self.assertEqual(payload["status"], "suppressed")
            self.assertFalse((stale / "closest_species_tree.svg").exists())
            self.assertTrue((stale / "phylogenetic_placement.json").is_file())

    def test_possible_carryover_no_longer_forces_mixed_signal_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database = self._database(root)
            species_root = root / "species_identification"
            self._sample(
                species_root,
                "salmonella_donor",
                "28901",
                "Salmonella enterica",
                [signal("28901", "Salmonella enterica", 100_000, 100.0)],
            )
            self._sample(
                species_root,
                "ecoli_recipient",
                "562",
                "Escherichia coli",
                [
                    signal("562", "Escherichia coli", 10_000, 95.0),
                    signal("28901", "Salmonella enterica", 500, 5.0),
                ],
            )
            results = reconcile_batch_identifications(
                species_output_root=species_root,
                database_root=database,
                expected_organisms={
                    "salmonella_donor": "Salmonella enterica",
                    "ecoli_recipient": "Escherichia coli",
                },
                platform="Nanopore",
            )
            recipient = results["ecoli_recipient"]
            self.assertEqual(recipient.status, "confirmed_species")
            self.assertEqual(recipient.reported_organism, "Escherichia coli")
            self.assertEqual(recipient.taxonomic_confirmation, "Pass")
            self.assertEqual(
                next(
                    item.status
                    for item in recipient.cross_sample_evidence
                    if item.species_taxid == "28901"
                ),
                "possible_cross_talk",
            )
            self.assertTrue(
                (species_root / "ecoli_recipient/cross_sample_signals.tsv").is_file()
            )
            provenance = __import__("json").loads(
                (species_root / "ecoli_recipient/analysis_provenance.json").read_text()
            )
            self.assertTrue(provenance["cohort_review"]["performed"])
            self.assertEqual(provenance["result_status"], "confirmed_species")

            summary = root / "temp_summary.tsv"
            summary.write_text(
                "Sample\tMean quality\tExpected organism\tIdentified organism\t"
                "% Match\ttax_confirm\n"
                "ecoli_recipient\t20\tEscherichia coli\tNo confident identification\t"
                "N/A\tFail\n"
            )
            update_summary_identifications(summary, results)
            with summary.open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["Identified organism"], "Escherichia coli")
            self.assertEqual(row["tax_confirm"], "Pass")

            final_summary = root / "final_summary.tsv"
            final_summary.write_text(
                "BactiPipe>>v0.1.0\n"
                "Sample ID\tIdentified organism\t% Match\tOverall Quality\n"
                "ecoli_recipient\tNo confident identification\tN/A\tFail\n"
            )
            update_summary_identifications(final_summary, results)
            with final_summary.open() as handle:
                lines = handle.read().splitlines()
            final_row = next(csv.DictReader(lines[1:], delimiter="\t"))
            self.assertEqual(final_row["Identified organism"], "Escherichia coli")
            self.assertNotIn("tax_confirm", final_row)

    def test_checkm_does_not_replace_independent_assembly_support(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database = self._database(root)
            species_root = root / "species_identification"
            self._sample(
                species_root,
                "suspect",
                "562",
                "Escherichia coli",
                [
                    signal("562", "Escherichia coli", 10_000, 70.0),
                    signal("28901", "Salmonella enterica", 4_000, 30.0),
                ],
            )

            review = reconcile_batch_identifications(
                species_output_root=species_root,
                database_root=database,
                expected_organisms={"suspect": "Escherichia coli"},
                platform="Illumina",
                checkm_contamination_by_sample={"suspect": 1.0},
            )["suspect"]
            self.assertEqual(review.status, "mixed_signal_review")

            contaminated = reconcile_batch_identifications(
                species_output_root=species_root,
                database_root=database,
                expected_organisms={"suspect": "Escherichia coli"},
                platform="Illumina",
                checkm_contamination_by_sample={"suspect": 99.84},
            )["suspect"]
            self.assertEqual(
                contaminated.status, "mixed_signal_review"
            )
            self.assertEqual(contaminated.reported_organism, "Escherichia coli")
            self.assertEqual(contaminated.match_percent, 99.0)

            # A genuinely corroborated pair remains Fail after reconciliation;
            # the independently computed evidence must survive JSON publication.
            import json
            from dataclasses import asdict
            from bactipipe.assembly_mixture import POLICY
            result_path = species_root / 'suspect' / 'species_identification.json'
            saved = json.loads(result_path.read_text())
            saved['mixed_signal']['assembly_evidence'] = {
                'policy': asdict(POLICY), 'status': 'assessed', 'primary_taxid': '562',
                'species': {'562': {'independent_support': True},
                            '28901': {'independent_support': True}}}
            result_path.write_text(json.dumps(saved))
            independent = reconcile_batch_identifications(
                species_output_root=species_root, database_root=database,
                expected_organisms={'suspect': 'Escherichia coli'}, platform='Illumina',
                checkm_contamination_by_sample={'suspect': 1.0})['suspect']
            self.assertEqual(independent.status, 'suspected_culture_contamination')
            self.assertEqual(independent.mixed_signal.assembly_evidence['policy']['version'], POLICY.version)

    def test_final_summary_verdict_is_recalculated_from_updated_temporary_row(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            temporary_summary = root / "temp.tsv"
            temporary_summary.write_text(
                "Sample\tMean quality\tqc_verdict\tExpected organism\t"
                "Identified organism\t% Match\tCoverage\tmin_cov\t"
                "cov_verdict\ttax_confirm\n"
                "S1\t20\tPass\tEscherichia coli\tEscherichia coli\t"
                "99.0\t100X\t30\tPass\tPass\n"
            )
            final_summary = root / "final.tsv"
            final_summary.write_text(
                "BactiPipe>>v0.1.0\n"
                "Sample ID\tMean quality\tExpected organism\tIdentified organism\t"
                "% Match\tCoverage Depth\tRequired depth\tCheckM completeness\t"
                "CheckM contamination\tOverall Quality\n"
                "S1\t20\tEscherichia coli\tEscherichia coli\t99.0\t100X\t"
                "30X\t99.0\t1.0\tFail\n"
            )
            update_final_summary_verdicts(final_summary, temporary_summary)
            _prefix, _fields, rows = __import__(
                "bactipipe.batch_species", fromlist=["_read_summary"]
            )._read_summary(final_summary)
            self.assertEqual(rows[0]["Overall Quality"], "Pass")

    @staticmethod
    def _database(root: Path) -> Path:
        database = root / "database"
        refs = database / "references"
        refs.mkdir(parents=True)
        rows = []
        for accession, taxid, name in (
            ("GCF_ECOLI", "562", "Escherichia coli"),
            ("GCF_SALM", "28901", "Salmonella enterica"),
        ):
            fasta = refs / f"{accession}.fna"
            fasta.write_text(f">{accession}\nACGTACGT\n")
            rows.append(
                {
                    "accession": accession,
                    "species_taxid": taxid,
                    "species_name": name,
                    "fasta_path": f"references/{accession}.fna",
                    "source": "RefSeq",
                    "reference_role": "type",
                    "assembly_level": "complete",
                    "sequence_origin": "isolate",
                    "type_material": "yes",
                    "canonical_for_tree": "yes",
                    "callable_for_identification": "yes",
                    "include_in_kraken": "yes",
                    "reporting_group": "",
                    "ncbi_taxonomy_name": name,
                    "gtdb_taxonomy": "",
                    "sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
                }
            )
        manifest = database / "reference_manifest.tsv"
        with manifest.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        return database

    @staticmethod
    def _sample(
        root: Path,
        sample_id: str,
        taxid: str,
        species_name: str,
        signals: list[TaxonomicSignal],
    ) -> None:
        sample = root / sample_id
        sample.mkdir(parents=True)
        (sample / "skani.tsv").write_text(
            "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\t"
            "Align_fraction_query\tRef_name\tQuery_name\n"
            f"GCF_{'ECOLI' if taxid == '562' else 'SALM'}.fna\tquery.fna\t"
            "99.0\t95.0\t95.0\tref\tquery\n"
        )
        write_taxonomic_signals(signals, sample / "taxonomic_signals.tsv")
        initial = decide_identification(
            sample_id=sample_id,
            expected_organism=species_name,
            ani_hits=(),
        )
        write_identification_outputs(initial, sample)


if __name__ == "__main__":
    unittest.main()
