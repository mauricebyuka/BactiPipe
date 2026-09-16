"""Bounded, isolated Illumina QC/assembly jobs; batch reporting stays in the parent."""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path


def available_cpus():
    count = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1)
    try:
        quota, period = Path('/sys/fs/cgroup/cpu.max').read_text().split()
        if quota != 'max':
            count = min(count, max(1, int(quota) // int(period)))
    except (OSError, ValueError):
        pass
    return max(1, count)


def cpu_budget(requested=None):
    if requested is not None and int(requested) < 1:
        raise ValueError('CPU budget must be positive')
    return min(int(requested or available_cpus()), available_cpus())


def available_memory_gb():
    import psutil
    # CPU-partitioned web workers share host RAM: reserve a proportional share
    # and headroom rather than letting each run claim all currently free memory.
    share = min(1.0, available_cpus() / max(1, os.cpu_count() or 1))
    usable = psutil.virtual_memory().available * share
    try:
        maximum = Path('/sys/fs/cgroup/memory.max').read_text().strip()
        if maximum != 'max':
            current = int(Path('/sys/fs/cgroup/memory.current').read_text())
            usable = min(usable, max(0, int(maximum) - current))
    except (OSError, ValueError):
        pass
    return max(0, usable * 0.8 / (1024 ** 3))


@dataclass(frozen=True)
class Resources:
    total_cpus: int
    jobs: int
    threads_per_sample: int
    memory_per_sample_gb: int
    memory_budget_gb: float


def plan_resources(samples, total_cpus, *, maximum_jobs=3, maximum_threads=20,
                   memory_gb=None, assembler='Skesa'):
    if min(samples, total_cpus, maximum_jobs, maximum_threads) < 1:
        raise ValueError('Samples, CPUs and concurrency limits must be positive')
    budget = available_memory_gb() if memory_gb is None else float(memory_gb)
    if not math.isfinite(budget) or budget < 4:
        raise ValueError('Less than 4 GiB available for Illumina assembly; wait for resources')
    # 24 GiB is a scheduling estimate, not a guarantee that every input fits.
    # On smaller hosts one job receives the available bounded budget.
    memory_jobs = max(1, int(budget // 24))
    cpu_jobs = max(1, total_cpus // min(8, maximum_threads))
    jobs = min(samples, maximum_jobs, cpu_jobs, memory_jobs)
    if assembler.lower() not in {'skesa', 'spades'}:
        jobs = 1  # Unicycler behavior is not part of this rollout.
        maximum_threads = total_cpus
    return Resources(total_cpus, jobs, min(maximum_threads, total_cpus // jobs),
                     max(1, min(32, int(budget // jobs))), round(budget, 2))


@dataclass(frozen=True)
class SampleJob:
    sample: str
    reads: tuple[str, ...]
    qc_root: str
    assembly_root: str
    assembler: str
    assembler_version: str
    minimum_quality: float
    threads: int
    memory_gb: int


def atomic_json(path, payload):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(payload, indent=2) + '\n')
    temporary.replace(path)


def prepare_sample(job):
    from .progress import emit
    from bactipipe.illumina_filtering import select_reads, assemble_selected
    root = Path(job.qc_root) / job.sample
    root.mkdir(parents=True, exist_ok=True)
    log = str(root / 'illumina_worker.log')
    state = root / 'illumina_worker_state.json'
    result = dict(sample=job.sample, selected=None, genome=None, failure=None)
    stage = 'qc'
    try:
        input_bytes = sum(Path(p).stat().st_size for p in job.reads)
    except OSError:
        input_bytes = None
    emit('read_qc', sample=job.sample, metrics={'input_bytes': input_bytes,
        'threads': job.threads, 'memory_limit_gib': job.memory_gb})
    try:
        atomic_json(state, dict(sample=job.sample, status='running', stage=stage))
        selected = select_reads(job.reads, job.qc_root, job.sample, job.minimum_quality, job.threads, log)
        result['selected'] = selected
        emit('read_qc', 'completed', sample=job.sample, metrics={
            'selected_bases': selected.get('selected', {}).get('total_bases'), 'selected_reads': selected.get('selected', {}).get('total_reads'),
            'raw_bases': selected.get('raw', {}).get('total_bases'), 'raw_reads': selected.get('raw', {}).get('total_reads'),
            'filter_applied': selected.get('filter_applied')})
        if selected['selected']['mean_quality'] >= job.minimum_quality:
            stage = 'assembly'
            emit('assembly', sample=job.sample, metrics={'threads': job.threads, 'memory_limit_gib': job.memory_gb})
            atomic_json(state, dict(sample=job.sample, status='running', stage=stage))
            result['genome'] = assemble_selected(
                job.sample, selected, job.assembly_root, job.assembler,
                job.assembler_version, job.threads, log,
                memory_gb=job.memory_gb if job.assembler.lower() in {'skesa', 'spades'} else None,
            )
            emit('assembly', 'completed' if result['genome'] else 'failed', sample=job.sample)
    except Exception as exc:
        emit('read_qc' if stage == 'qc' else 'assembly', 'failed', sample=job.sample)
        result['failure'] = dict(stage=stage, reason=f'{type(exc).__name__}: {exc}')
    atomic_json(state, dict(sample=job.sample, status='failed' if result['failure'] else 'completed',
                           stage=stage, failure=result['failure']))
    return result


def _launch(job):
    root = Path(job.qc_root) / job.sample
    root.mkdir(parents=True, exist_ok=True)
    specification, result_path = root / 'illumina_job.json', root / 'illumina_preparation.json'
    atomic_json(specification, asdict(job))
    result_path.unlink(missing_ok=True)  # Never accept a stale successful worker result.
    # These caps are for Python QC; assemble() overrides OpenMP only in its
    # assembler child to match the job's allocated --threads/--cores budget.
    env = dict(os.environ, MPLBACKEND='Agg', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    # Separate interpreters avoid the Python GIL and matplotlib thread unsafety.
    # Inherit the pipeline process group so Stop analysis reaches all children.
    with (root / 'illumina_worker.stdout.log').open('w') as log:
        from .progress import PREFIX
        proc = subprocess.Popen([sys.executable, '-m', 'bactipipe.illumina_parallel', str(specification)],
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env, text=True,
                                errors='replace', bufsize=1)
        for line in proc.stdout:
            log.write(line)
            if line.startswith(PREFIX):
                try:
                    os.write(1, line.encode())
                except OSError:
                    pass
        proc.wait()
        if hasattr(proc.stdout, 'close'):
            proc.stdout.close()
    if proc.returncode or not result_path.is_file():
        raise RuntimeError(f'QC/assembly worker exited {proc.returncode}; see {root / "illumina_worker.stdout.log"}')
    result = json.loads(result_path.read_text())
    if result.get('sample') != job.sample:
        raise RuntimeError('QC/assembly worker returned an unexpected sample')
    return result


def prepare_batch(jobs, resources, report=lambda message: None, *, launch=None):
    from .progress import emit
    jobs = tuple(jobs)
    emit('assembly', total=len(jobs), metrics={'concurrent_samples': resources.jobs,
         'threads_per_sample': resources.threads_per_sample, 'memory_per_sample_gib': resources.memory_per_sample_gb,
         'allocated_cpus': resources.total_cpus, 'memory_budget_gib': resources.memory_budget_gb})
    if len({job.sample for job in jobs}) != len(jobs):
        raise ValueError('Duplicate sample IDs cannot be analyzed concurrently')
    launch = launch or _launch
    outcomes = {}
    report(f'Illumina preparation: {resources.jobs} concurrent samples, '
           f'{resources.threads_per_sample} CPUs and {resources.memory_per_sample_gb} GiB per sample')
    report(f'Illumina QC/assembly allocation: {resources.jobs} concurrent samples x '
           f'{resources.threads_per_sample} threads = {resources.jobs * resources.threads_per_sample}/'
           f'{resources.total_cpus} allocated CPUs; {len(jobs)} samples total. '
           'Read QC/filtering and assembly are included; identification follows preparation.')
    with ThreadPoolExecutor(max_workers=resources.jobs) as executor:
        pending = {executor.submit(launch, job): job for job in jobs}
        for future in as_completed(pending):
            job = pending[future]
            try:
                outcomes[job.sample] = future.result()
            except Exception as exc:
                outcomes[job.sample] = dict(sample=job.sample, selected=None, genome=None,
                    failure=dict(stage='preparation', reason=f'{type(exc).__name__}: {exc}'))
            failure = outcomes[job.sample]['failure']
            emit('read_qc', 'completed' if outcomes[job.sample].get('selected') else 'failed', sample=job.sample)
            emit('assembly', 'completed' if outcomes[job.sample].get('genome') else ('failed' if failure and failure.get('stage') == 'assembly' else 'skipped'), sample=job.sample)
            report(f"Illumina preparation progress: {len(outcomes)}/{len(jobs)} finished; "
                   f"sample {job.sample}: QC/assembly {'failed: ' + failure['reason'] if failure else 'completed'}")
    # Parent consumes outcomes in sample-sheet order, not completion order.
    return {job.sample: outcomes[job.sample] for job in jobs}


if __name__ == '__main__':
    specification = Path(sys.argv[1])
    payload = json.loads(specification.read_text())
    payload['reads'] = tuple(payload['reads'])
    atomic_json(specification.with_name('illumina_preparation.json'), prepare_sample(SampleJob(**payload)))
