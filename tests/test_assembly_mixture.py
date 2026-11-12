import json
import tempfile
import unittest
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import patch

from bactipipe.assembly_mixture import (
    POLICY, assess_paf, apply_assembly_evidence, independent_mixture_supported,
    interval_bases, investigate_assembly,
)
from bactipipe.species_identification import (
    AniHit, TaxonomicSignal, assess_mixed_signals, decide_identification,
    checkm_contamination_evidence, ReferenceCatalog, ReferenceRecord,
)
from bactipipe.species_workflow import legacy_summary_fields
from bactipipe.qc_status import interpret


def mixed():
    return assess_mixed_signals([
        TaxonomicSignal("1", "Primary", 10000, 65, 95, 10000),
        TaxonomicSignal("2", "Secondary", 5000, 32, 70, 5000)])


def evidence(species):
    return {"policy": asdict(POLICY), "status": "assessed",
            "primary_taxid": "1", "species": species,
            "reference_comparisons": [{"secondary_taxid": "2", "measurements": [
                {"ANI": "99.08", "Align_fraction_ref": "89.65", "Align_fraction_query": "89.13"}]}]}


def supported():
    return evidence({"1": {"independent_support": True},
                     "2": {"independent_support": True}})


def paf(query, target, position, matches=9900, qlen=10000, end=None):
    end = qlen if end is None else end
    return f"{query}\t{qlen}\t0\t{end}\t+\t{target}\t1000000\t{position}\t{position+end}\t{matches}\t{end}\t60\n"


class AssemblyMixtureTests(unittest.TestCase):
    def measure(self, lines):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "alignments.paf"
            path.write_text("".join(lines))
            return assess_paf(path, targets={"a": {"taxid": "1", "length": 1000000},
                "b": {"taxid": "2", "length": 1000000}}, primary_taxid="1")

    def test_shared_windows_cannot_create_two_genomes(self):
        data = self.measure([line for i in range(30) for line in
            [paf(f"q{i}", "a", i * 30000), paf(f"q{i}", "b", i * 30000, 9870)]])
        self.assertEqual(data["2"]["exclusive_bases"], 0)
        self.assertEqual(data["2"]["shared_with_primary_percent"], 100)
        result = apply_assembly_evidence(mixed(), evidence(data))
        self.assertEqual(result.status, "reference_ambiguity")
        self.assertFalse(independent_mixture_supported(result))

    def test_distinct_distributed_windows_support_mixture_on_both_platforms(self):
        data = self.measure([line for i in range(20) for line in
            [paf(f"a{i}", "a", i * 45000), paf(f"b{i}", "b", i * 45000)]])
        result = apply_assembly_evidence(mixed(), evidence(data))
        self.assertTrue(independent_mixture_supported(result))
        for platform in ("Illumina", "Nanopore"):
            decision = self.decide(result, platform=platform)
            self.assertEqual(decision.status, "suspected_culture_contamination")
            self.assertEqual(decision.reported_organism, "Primary")
            self.assertIn("do not establish its source", decision.reason)

    def test_repeated_same_locus_is_not_independent_support(self):
        data = self.measure([paf(f"q{i}", "b", 1000) for i in range(200)])
        self.assertEqual(data["2"]["exclusive_bases"], 10000)
        self.assertFalse(data["2"]["independent_support"])

    def test_single_accessory_region_does_not_pass_distribution(self):
        data = self.measure([paf(f"q{i}", "b", i * 10000) for i in range(10)])
        self.assertEqual(data["2"]["exclusive_bases"], 100000)
        self.assertFalse(data["2"]["independent_support"])

    def test_partial_and_low_identity_alignments_do_not_support_genome(self):
        data = self.measure([line for i in range(20) for line in
            [paf(f"a{i}", "a", i * 45000, matches=9400),
             paf(f"b{i}", "b", i * 45000, matches=4900, end=5000)]])
        self.assertFalse(any(row["independent_support"] for row in data.values()))

    def test_tied_reference_evidence_does_not_change_species_names(self):
        data = self.measure([line for i in range(20) for line in
            [paf(f"q{i}", "a", i * 40000), paf(f"q{i}", "b", i * 40000)]])
        decision = self.decide(apply_assembly_evidence(mixed(), evidence(data)))
        self.assertEqual(decision.status, "confirmed_species")
        self.assertEqual(decision.reported_organism, "Primary")
        self.assertEqual(decision.mixed_signal.status, "reference_ambiguity")

    def test_missing_or_old_evidence_never_proves_contamination(self):
        for platform in ("Illumina", "Nanopore"):
            for checkm in (0.37, 99):
                result = self.decide(mixed(), platform=platform, checkm=checkm)
                self.assertEqual(result.status, "mixed_signal_review")
                self.assertEqual(result.reported_organism, "Primary")
                fields = legacy_summary_fields(result)
                self.assertEqual(fields.identified_organism, "Primary")
                self.assertEqual(fields.taxonomic_confirmation, "MixedSignalReview")
        old = supported(); old["policy"]["version"] = "old"
        self.assertFalse(independent_mixture_supported(apply_assembly_evidence(mixed(), old)))

    def test_cross_talk_exclusion_cannot_be_revived_by_saved_evidence(self):
        result = apply_assembly_evidence(replace(mixed(), status="single_signal", qualifying_secondary=()), supported())
        self.assertFalse(independent_mixture_supported(result))
        self.assertEqual(self.decide(result).status, "confirmed_species")

    def test_changed_primary_does_not_reuse_comparisons(self):
        result = apply_assembly_evidence(replace(mixed(), primary=TaxonomicSignal("3", "New", 10000, 70)), supported())
        self.assertFalse(independent_mixture_supported(result))

    def test_serialization_retains_policy_and_evidence(self):
        result = self.decide(apply_assembly_evidence(mixed(), supported()))
        saved = json.loads(json.dumps(result.to_dict()))
        self.assertEqual(saved["mixed_signal"]["assembly_evidence"]["policy"]["version"], POLICY.version)

    def test_shared_regions_without_reference_comparison_stay_review(self):
        data = self.measure([line for i in range(20) for line in
            [paf(f"q{i}", "a", i * 40000), paf(f"q{i}", "b", i * 40000)]])
        saved = evidence(data); saved["reference_comparisons"] = []
        self.assertEqual(self.decide(apply_assembly_evidence(mixed(), saved)).status, "mixed_signal_review")

    def test_review_does_not_erase_measurement_failures(self):
        self.assertEqual(interpret("Unknown", "Primary", ["Pass"] * 4, isolate_qc_review=True)["overall"], "Review")
        self.assertEqual(interpret("Unknown", "Primary", ["Pass", "Fail", "Pass", "Pass"], isolate_qc_review=True)["overall"], "Fail")
        self.assertEqual(interpret("Unknown", "Primary", ["Pass", "Pass", "Pass", "Fail"], isolate_qc_review=True)["overall"], "Fail")

    def test_interval_union(self):
        self.assertEqual(interval_bases([(0, 100), (50, 150), (20, 30), (200, 300)]), 250)

    def test_reference_or_tool_failure_is_unassessed_not_clean(self):
        for header, size, expect_command in (("plasmid", 600000, False),
                                             ("chromosome", 10000, False),
                                             ("chromosome", 600000, True)):
            with self.subTest(header=header, size=size), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                assembly = root / 'assembly.fasta'; assembly.write_text('>a\n' + 'ACGT' * 3000 + '\n')
                records = []
                for taxid in ('1', '2'):
                    path = root / f'{taxid}.fna'
                    path.write_text(f'>{header}\n' + 'A' * size + '\n')
                    records.append(ReferenceRecord(taxid, taxid, taxid, path.name,
                        'test', 'representative', 'complete', 'isolate', False, True, True))
                catalog = ReferenceCatalog(records, root=root)
                with patch('bactipipe.assembly_mixture.subprocess.run', side_effect=FileNotFoundError('tool unavailable')) as run:
                    data = investigate_assembly(assembly=assembly, catalog=catalog,
                        mixed=mixed(), output_directory=root / 'output')
                self.assertEqual(run.called, expect_command)
                self.assertEqual(data['status'], 'not_assessed')
                self.assertEqual(self.decide(apply_assembly_evidence(mixed(), data)).status, 'mixed_signal_review')
                self.assertTrue((root / 'output/assembly_mixture.json').is_file())
                self.assertFalse((root / 'output/candidate_chromosomes.fna').exists())
                self.assertFalse((root / 'output/assembly_windows.fna').exists())

    @staticmethod
    def decide(assessment, platform="Nanopore", checkm=0.37):
        return decide_identification(sample_id="S", platform=platform,
            ani_hits=[AniHit("ref1", "1", "Primary", "r.fna", 99.45, 95.21, 90.53)],
            mixed_signal=assessment, cohort_review_complete=True,
            contamination_evidence=checkm_contamination_evidence(checkm))
