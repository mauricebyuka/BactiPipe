"""Bounded Nanopore assembly allocation and completion-order progress."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import math
from pathlib import Path

from .illumina_parallel import available_memory_gb, atomic_json
from .progress import emit


@dataclass(frozen=True)
class AssemblyResources:
    total_cpus: int
    jobs: int
    threads_per_sample: int
    memory_budget_gb: float
    estimated_memory_per_sample_gb: int = 8


def plan_flye_resources(samples, total_cpus, *, requested_threads=None, memory_gb=None):
    """At most six jobs and 10 default threads; explicit thread overrides survive.

    The existing Nanopore 8 GiB/job estimate is not an imposed Flye memory limit.
    CPU affinity/cgroup and shared-host memory limits are resolved by the caller.
    """
    if samples < 1 or total_cpus < 1 or (requested_threads is not None and requested_threads < 1):
        raise ValueError('Sample and CPU allocations must be positive')
    memory = available_memory_gb() if memory_gb is None else float(memory_gb)
    if not math.isfinite(memory) or memory < 4:
        raise ValueError('Less than 4 GiB available for assembly; wait for resources')
    threads = min(total_cpus, requested_threads or 10)
    jobs = min(samples, 6, max(1, int(memory // 8)), max(1, total_cpus // threads))
    return AssemblyResources(total_cpus, jobs, threads, round(memory, 2))


def run_assembly_jobs(jobs, *, workers, worker, report, state_path, executor_factory=ProcessPoolExecutor):
    """A slow early sample must not hide completions from later samples."""
    jobs = tuple(jobs)
    if workers < 1 or len({job.sample for job in jobs}) != len(jobs):
        raise ValueError('Require positive worker count and unique samples')
    state_path = Path(state_path)
    state_path.parent.mkdir(parents=True, exist_ok=True)
    state = {'stage': 'assembly', 'total': len(jobs), 'finished': 0, 'completed': 0,
             'failed': 0, 'samples': {job.sample: 'unfinished' for job in jobs}}
    atomic_json(state_path, state)
    with executor_factory(max_workers=workers) as executor:
        futures = {executor.submit(worker, job): job for job in jobs}
        for future in as_completed(futures):
            job = futures[future]
            reason = ''
            try:
                future.result()
            except Exception as exc:
                reason = f'Assembly worker failed: {exc}'
            genome = Path(job.assembly_dir) / 'genomes' / f'{job.sample}.fasta'
            complete = not reason and genome.is_file() and genome.stat().st_size > 0
            if not complete:
                failure = Path(job.assembly_dir).parent / 'qc_out' / job.sample / 'assembly_failure.txt'
                failure.parent.mkdir(parents=True, exist_ok=True)
                if reason or not failure.exists():
                    failure.write_text((reason or 'Assembly did not produce a usable FASTA.') + '\n')
            state['finished'] += 1
            state['completed' if complete else 'failed'] += 1
            state['samples'][job.sample] = 'completed' if complete else 'failed'
            emit('assembly', 'completed' if complete else 'failed', sample=job.sample)
            atomic_json(state_path, state)
            report(f"Assembly progress: {state['finished']}/{len(jobs)} finished "
                   f"({state['completed']} completed, {state['failed']} failed); "
                   f"sample {job.sample}: {'completed' if complete else 'failed'}")
    return state
