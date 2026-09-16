import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from bactipipe.adaptive_assembly import (
    AssemblyPolicy, adaptive_flye, fastq_bases, retain_assembly, retry_decision,
)


class AdaptiveAssemblyTests(unittest.TestCase):
    def decide(self, **kwargs):
        values = dict(success=False, exit_code=1, log_text="", selected_bases=250_000_000,
                      assembly_bases=0, more_reads=True)
        values.update(kwargs)
        return retry_decision(**values)

    def test_low_coverage_failure_only(self):
        self.assertTrue(self.decide(log_text="Overlap-based coverage: 12\nNo disjointigs were assembled")[0])
        self.assertFalse(self.decide(log_text="Overlap-based coverage: 1173\nNo disjointigs were assembled")[0])
        self.assertFalse(self.decide(log_text="No disjointigs were assembled")[0])
        self.assertFalse(self.decide(exit_code=124, log_text="Overlap-based coverage: 12\nNo disjointigs")[0])
        self.assertFalse(self.decide(log_text="Overlap-based coverage: 12\nNo disjointigs", more_reads=False)[0])

    def test_successful_sampling(self):
        self.assertFalse(self.decide(success=True, assembly_bases=5_000_000)[0])
        self.assertTrue(self.decide(success=True, assembly_bases=8_000_000)[0])
        self.assertFalse(self.decide(success=True, assembly_bases=8_000_000, more_reads=False)[0])

    def test_defaults(self):
        self.assertEqual(AssemblyPolicy().budgets, (250_000_000, 500_000_000, 1_000_000_000))

    def test_fastq_and_fasta_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.fastq"
            path.write_text("@x\nACGT\n+\nIIII\n")
            self.assertEqual(fastq_bases(path), 4)
            path.write_text("@x\nACGT\n+\nIII\n")
            with self.assertRaises(ValueError):
                fastq_bases(path)
            path.write_text(">one\nACGT\n>two\nNNN\n>short\nA\n")
            result = retain_assembly(path, Path(tmp) / "assembly.fasta", 3)
            self.assertEqual(result["assembly_bases"], 7)
            self.assertEqual(result["contigs"], 2)

    def run_fixture(self, root, *, fail_second=False):
        source = root / "reads.fastq"
        source.write_text("@x\n" + "A" * 100 + "\n+\n" + "I" * 100 + "\n")
        self.commands = []

        def fake(argv, log, *, output=None, timeout=None):
            self.commands.append(argv)
            Path(log).write_text("")
            code = 0
            if argv[0] == "filtlong":
                count = int(argv[argv.index("--target_bases") + 1])
                Path(output).write_text("@x\n" + "A" * count + "\n+\n" + "I" * count + "\n")
            else:
                out = Path(argv[argv.index("--out-dir") + 1])
                out.mkdir(exist_ok=True)
                if fail_second and "attempt-02" in str(out):
                    code = 1
                    Path(log).write_text("tool failed")
                else:
                    (out / "assembly.fasta").write_text(">one\nAAAA\n")
            return {"argv": argv, "exit_code": code, "seconds": .01}

        with patch("bactipipe.adaptive_assembly.subprocess.check_output", return_value="v1"), \
                patch("bactipipe.adaptive_assembly.run_tool", side_effect=fake):
            result = adaptive_flye(reads=source, output=root / "final.fasta", audit_dir=root / "audit",
                                   cpus=2, policy=AssemblyPolicy(budgets=(20, 40, 80), sampling_depth=10, minimum_contig=1))
            calls = len(self.commands)
            resumed = adaptive_flye(reads=source, output=root / "final.fasta", audit_dir=root / "audit",
                                    cpus=2, policy=AssemblyPolicy(budgets=(20, 40, 80), sampling_depth=10, minimum_contig=1))
            self.assertEqual(calls, len(self.commands))
            self.assertEqual(result, resumed)
        self.assertEqual(fastq_bases(source), 100)
        self.assertEqual(list((root / "audit").rglob("selected.fastq")), [])
        for cmd in self.commands:
            if cmd[0] == "flye":
                self.assertNotIn("--g", cmd)
                self.assertNotIn("--asm-coverage", cmd)
        return result

    def test_adaptive_retry_and_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = self.run_fixture(Path(tmp))
            self.assertEqual(result["selected_attempt"], 2)
            self.assertEqual(len(result["attempts"]), 2)

    def test_preserve_success_after_later_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = self.run_fixture(Path(tmp), fail_second=True)
            self.assertEqual(result["selected_attempt"], 1)
            self.assertEqual(len(result["attempts"]), 2)

    def test_cli_passes_actual_length_option(self):
        path = Path(__file__).parents[1] / "bactipipe/scripts/run_nanopore_qc.py"
        import ast
        tree = ast.parse(path.read_text())
        jobs = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name) and node.func.id == "NanoporeAssemblyJob"]
        self.assertEqual(len(jobs), 1)
        args = {kw.arg: ast.unparse(kw.value) for kw in jobs[0].keywords}
        self.assertEqual(args["min_length"], "args.min_read_length")


if __name__ == "__main__":
    unittest.main()
