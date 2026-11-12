"""Bounded per-sample nomination/ANI; reporting remains ordered in the parent.

Reference indexes are read-only. Each callback owns a separate sample workspace.
The memory estimate is a scheduling guard, not a measured peak or hard RAM limit.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, asdict
import json
import math
from pathlib import Path

from .illumina_parallel import available_memory_gb, cpu_budget
from .progress import emit


@dataclass(frozen=True)
class IdentificationResources:
    jobs: int
    threads_per_sample: int
    memory_budget_gb: float
    estimated_memory_per_sample_gb: float


def plan_identification(samples, cpus, database_root, *, memory_gb=None):
    budget = available_memory_gb() if memory_gb is None else float(memory_gb)
    if not math.isfinite(budget) or budget < 0:
        raise ValueError('Invalid identification memory budget')
    # Account conservatively for separate KMA processes, even when pages can be
    # shared by the OS. Unknown/legacy index layouts stay serial.
    estimate, known_index = 32.0, False
    try:
        identity = json.loads((Path(database_root) / 'database_identity.json').read_text())
        components = identity.get('components', {})
        sizes = [int(item['size_bytes']) for item in components.values()]
        known_index = identity.get('database_kind') == 'speciesfinder_kma' and bool(sizes) and min(sizes) >= 0
        if known_index:
            estimate = max(32.0, 2 * sum(sizes) / (1024 ** 3) + 4)
    except (OSError, ValueError, TypeError, KeyError, AttributeError):
        pass
    cpus = cpu_budget(cpus)
    jobs = min(max(1, samples), 2 if known_index else 1,
               max(1, cpus // 8), max(1, int(budget // estimate)))
    return IdentificationResources(jobs, max(1, cpus // jobs), round(budget, 2), round(estimate, 2))


def identify_batch(jobs, identify, *, cpus, database_root, report, audit_path,
                   memory_gb=None):
    """jobs: sample -> (expected organism, assembly, reads). Return value/errors.

    Only expensive per-sample work runs concurrently. Batch reconciliation,
    summary writing, CheckM and publication remain in their existing order.
    """
    if not jobs:
        return {}
    resources = plan_identification(len(jobs), cpus, database_root, memory_gb=memory_gb)
    Path(audit_path).write_text(json.dumps(asdict(resources), indent=2) + '\n')
    report(f'Identification: {len(jobs)} samples; up to {resources.jobs} concurrent; '
           f'{resources.threads_per_sample} threads/sample; '
           f'{resources.memory_budget_gb:g} GiB memory budget '
           f'(estimated {resources.estimated_memory_per_sample_gb:g} GiB/sample).')

    def execute(sample, inputs):
        emit('identification', sample=sample, metrics={'threads': resources.threads_per_sample})
        report(f'Sample {sample}: identification started.')
        try:
            result = identify(sample, *inputs, resources.threads_per_sample)
        except Exception as exc:
            emit('identification', 'failed', sample=sample)
            report(f'Sample {sample}: identification failed ({type(exc).__name__}).')
            return exc
        emit('identification', 'completed', sample=sample)
        report(f'Sample {sample}: identification completed.')
        return result

    results = {}
    with ThreadPoolExecutor(max_workers=resources.jobs) as pool:
        futures = {pool.submit(execute, sample, inputs): sample for sample, inputs in jobs.items()}
        for future in as_completed(futures):
            results[futures[future]] = future.result()
    return results


def result_for(results, sample):
    result = results[sample]
    if isinstance(result, Exception):
        raise result
    return result
