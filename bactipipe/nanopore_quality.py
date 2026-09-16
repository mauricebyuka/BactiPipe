"""Probability-averaged Nanopore read quality and Filtlong's required units.

Q_read = -10 log10(mean(10 ** (-Q_base / 10))). This is not the
arithmetic mean of base Phred scores. Filtlong without an external reference
uses 100 * (1 - mean_error_probability) for its --min_mean_q threshold.
"""
import math


def filtlong_min_mean_identity(read_quality: float) -> float:
    q = float(read_quality)
    if not math.isfinite(q) or not 0 <= q <= 60:
        raise ValueError('Nanopore read quality must be finite and between 0 and 60')
    return 100.0 * (1.0 - 10.0 ** (-q / 10.0))


def mean_read_quality(base_qualities) -> float:
    qualities = list(base_qualities)
    if not qualities:
        raise ValueError('A read must contain quality scores')
    probability = sum(10.0 ** (-float(q) / 10.0) for q in qualities) / len(qualities)
    return -10.0 * math.log10(probability)
