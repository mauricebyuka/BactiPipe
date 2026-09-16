"""Bounded, size-independent Flye input selection for unknown Nanopore isolates.

The retry depth is a computational sampling target, NOT a diagnostic QC cutoff.
Full retained reads remain untouched; CheckM and final yield depth run afterwards.
"""
from dataclasses import asdict, dataclass
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

from bactipipe.nanopore_quality import filtlong_min_mean_identity


@dataclass(frozen=True)
class AssemblyPolicy:
    budgets: tuple = (250_000_000, 500_000_000, 1_000_000_000)
    sampling_depth: float = 50.0
    attempt_timeout: int = 3600
    minimum_contig: int = 400

    def __post_init__(self):
        if not self.budgets or len(self.budgets) > 3 or any(b <= 0 for b in self.budgets):
            raise ValueError("Use one to three positive assembly budgets")
        if tuple(sorted(set(self.budgets))) != tuple(self.budgets):
            raise ValueError("Assembly budgets must be strictly increasing")
        if self.sampling_depth <= 0 or self.attempt_timeout <= 0 or self.minimum_contig < 1:
            raise ValueError("Assembly policy limits must be positive")


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fastq_bases(path):
    """Validate four-line FASTQ without loading reads into memory."""
    opener = gzip.open if str(path).endswith(".gz") else open
    bases = 0
    with opener(path, "rb") as handle:
        while header := handle.readline():
            seq, plus, qual = (handle.readline().rstrip(b"\r\n") for _ in range(3))
            if not header.startswith(b"@") or not plus.startswith(b"+") or not seq or len(seq) != len(qual):
                raise ValueError("Malformed four-line FASTQ")
            bases += len(seq)
    if not bases:
        raise ValueError("No reads remain for assembly")
    return bases


def run_tool(argv, log, *, output=None, timeout=3600):
    """Keep pipeline process-group membership so Stop analysis still works."""
    import psutil
    started = time.monotonic()
    with Path(log).open("ab") as err:
        destination = Path(output).open("wb") if output else err
        try:
            process = subprocess.Popen(argv, stdout=destination, stderr=err)
            try:
                code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                # Kill only this command's descendants, never sibling sample jobs.
                try:
                    parent = psutil.Process(process.pid)
                    descendants = parent.children(recursive=True)
                    parent.suspend()
                    for child in descendants:
                        try:
                            child.kill()
                        except psutil.NoSuchProcess:
                            pass
                    parent.kill()
                    psutil.wait_procs(descendants, timeout=10)
                except psutil.NoSuchProcess:
                    pass
                process.wait()
                code = 124
        finally:
            if output:
                destination.close()
    return {"argv": argv, "exit_code": code, "seconds": round(time.monotonic() - started, 2)}


def retry_decision(*, success, exit_code, log_text, selected_bases, assembly_bases,
                   more_reads, sampling_depth=50.0):
    """Additional reads help under-sampling, not every assembly failure."""
    if success:
        if not more_reads:
            return False, "available_reads_exhausted"
        if selected_bases / assembly_bases < sampling_depth:
            return True, "assembly_subset_below_sampling_target"
        return False, "sampling_target_met"
    if exit_code in (124, 137, -9) or re.search(
        r"out of memory|cannot allocate|no space left|permission denied|not found", log_text, re.I
    ):
        return False, "tool_resource_or_timeout_failure"
    coverage = re.findall(r"Overlap-based coverage:\s*([\d.]+)", log_text, re.I)
    low_coverage = bool(coverage and 0 < float(coverage[-1]) < sampling_depth)
    assembly_failure = bool(re.search(r"no disjointigs|not enough (?:reads|overlaps)|insufficient.*coverage", log_text, re.I))
    if more_reads and low_coverage and assembly_failure:
        return True, "flye_low_overlap_coverage"
    return False, "assembly_failure_not_explained_by_low_coverage"


def retain_assembly(source, target, minimum_contig):
    """Apply the same contig cutoff as the existing Flye publication path."""
    records = []
    for line in Path(source).read_text().splitlines():
        if line.startswith(">"):
            records.append([line, ""])
        elif line.strip():
            if not records:
                raise ValueError("Invalid Flye FASTA")
            records[-1][1] += line.strip()
    records = [(head, seq) for head, seq in records if len(seq) >= minimum_contig]
    if not records:
        raise ValueError("Flye produced no reportable contigs")
    Path(target).write_text("".join(f"{head}\n{seq}\n" for head, seq in records))
    return {"assembly_bases": sum(len(seq) for _, seq in records), "contigs": len(records),
            "assembly_sha256": sha256(target)}


def adaptive_flye(*, reads, output, audit_dir, cpus, min_quality=15, min_length=1000,
                  policy=AssemblyPolicy(), report=lambda message: None):
    """Return an audit; atomically publish a valid assembly, or raise on failure.

    Checkpoints are bound to input bytes, code, tool versions and sampling policy.
    A later failed attempt never destroys a previously successful assembly.
    """
    reads, output, root = Path(reads).resolve(), Path(output), Path(audit_dir)
    if cpus < 1 or min_length < 1 or output.resolve() == reads:
        raise ValueError("Invalid adaptive assembly resources or output path")
    root.mkdir(parents=True, exist_ok=True)
    versions = {tool: subprocess.check_output([tool, "--version"], stderr=subprocess.STDOUT,
                text=True, timeout=30).strip() for tool in ("filtlong", "flye")}
    signature = {"input_sha256": sha256(reads), "policy": asdict(policy), "cpus": cpus,
                 "min_quality": min_quality, "min_length": min_length, "versions": versions,
                 "implementation_sha256": sha256(__file__)}
    signature = json.loads(json.dumps(signature))
    manifest = root / "manifest.json"
    if manifest.exists() and json.loads(manifest.read_text()) != signature:
        raise ValueError("Adaptive assembly inputs changed; use a fresh sample workspace")
    atomic_json(manifest, signature)
    available_bases = fastq_bases(reads)
    attempts = []
    best = None
    for index, budget in enumerate(policy.budgets):
        work = root / f"attempt-{index + 1:02d}"
        work.mkdir(exist_ok=True)
        record_path, fasta = work / "result.json", work / "assembly.fasta"
        if record_path.exists():
            record = json.loads(record_path.read_text())
            if record["success"] and (not fasta.exists() or sha256(fasta) != record["assembly_sha256"]):
                raise ValueError("Adaptive assembly checkpoint failed checksum verification")
        else:
            report(f"Unknown-organism Flye: selecting up to {budget:,} bases (attempt {index + 1})")
            selected = work / "selected.fastq"
            filt = run_tool(["filtlong", "--min_mean_q", str(filtlong_min_mean_identity(min_quality)),
                             "--min_length", str(min_length), "--target_bases", str(budget), str(reads)],
                            work / "filter.log", output=selected, timeout=policy.attempt_timeout)
            atomic_json(work / "filter_execution.json", filt)
            if filt["exit_code"]:
                selected.unlink(missing_ok=True)
                raise RuntimeError(f"Assembly read selection failed (exit {filt['exit_code']})")
            selected_bases = fastq_bases(selected)
            selected_hash = sha256(selected)
            # Filtlong can overshoot the budget by one read; identical selections add no evidence.
            if attempts and selected_hash == attempts[-1]["selected_sha256"]:
                selected.unlink()
                attempts[-1]["retry"] = False
                attempts[-1]["decision"] = "no_additional_eligible_reads"
                break
            flye_dir = work / "flye"
            argv = ["flye", "--nano-hq", str(selected), "--out-dir", str(flye_dir),
                    "--threads", str(cpus), "--iterations", "2"]
            if (flye_dir / "params.json").exists():
                argv.append("--resume")
            execution = run_tool(argv, work / "command.log", timeout=policy.attempt_timeout)
            record = {"budget_bases": budget, "selected_bases": selected_bases,
                      "selected_sha256": selected_hash, "execution": execution, "success": False,
                      "assembly_bases": 0, "attempt": index + 1}
            if execution["exit_code"] == 0:
                try:
                    record.update(retain_assembly(flye_dir / "assembly.fasta", fasta, policy.minimum_contig))
                    record["success"] = True
                except (OSError, ValueError) as error:
                    record["assembly_error"] = str(error)
            log_text = (work / "command.log").read_text(errors="replace")[-20000:]
            for name in ("flye.log", "assembly_info.txt"):
                if (flye_dir / name).is_file():
                    shutil.copyfile(flye_dir / name, work / name)
            if (work / "flye.log").is_file():
                log_text += (work / "flye.log").read_text(errors="replace")[-20000:]
            record["retry"], record["decision"] = retry_decision(
                success=record["success"], exit_code=execution["exit_code"], log_text=log_text,
                selected_bases=selected_bases, assembly_bases=record["assembly_bases"],
                more_reads=selected_bases < available_bases and index + 1 < len(policy.budgets),
                sampling_depth=policy.sampling_depth)
            atomic_json(record_path, record)
            # Only command-owned scratch in this exact attempt is removed. Keep all compact evidence.
            selected.unlink(missing_ok=True)
            if flye_dir.is_dir():
                shutil.rmtree(flye_dir)
        attempts.append(record)
        if record["success"]:
            # No completeness-based guessing: preserve earlier output if a retry loses >5% of sequence.
            if best is None or record["assembly_bases"] >= best["assembly_bases"] * .95:
                best = record
        report(f"Unknown-organism Flye attempt {index + 1}: {record['decision']}")
        if not record["retry"]:
            break
    result = {"schema_version": 1, "method": "adaptive_unknown_flye", "available_qc_bases": available_bases,
              "attempts": attempts, "selected_attempt": best["attempt"] if best else None,
              "status": "assembled" if best else "failed", "diagnostic_qc": "assessed_separately",
              "selection_policy": "latest_success_without_more_than_5_percent_sequence_loss"}
    atomic_json(root / "summary.json", result)
    if best is None:
        raise RuntimeError(f"Unknown-organism assembly failed: {attempts[-1]['decision']}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    shutil.copyfile(root / f"attempt-{best['attempt']:02d}" / "assembly.fasta", temporary)
    temporary.replace(output)
    return result
