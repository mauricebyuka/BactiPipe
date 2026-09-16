import math
import pytest
from bactipipe.nanopore_quality import filtlong_min_mean_identity, mean_read_quality


@pytest.mark.parametrize('q', [0, 7, 10, 15, 20, 28, 40, 60])
def test_filtlong_conversion_preserves_probability_based_q(q):
    identity = filtlong_min_mean_identity(q)
    assert -10 * math.log10(1 - identity / 100) == pytest.approx(q)


def test_nanopore_quality_is_not_arithmetic_base_quality():
    qualities = [2, 40, 40, 40]
    assert sum(qualities)/len(qualities) > 30
    assert mean_read_quality(qualities) < 9
    assert 100 * (1 - sum(10**(-q/10) for q in qualities)/4) < filtlong_min_mean_identity(15)
    assert mean_read_quality([20]*1000) == pytest.approx(20)


@pytest.mark.parametrize('value', [-1, 61, float('nan'), float('inf')])
def test_invalid_quality_is_rejected(value):
    with pytest.raises(ValueError):
        filtlong_min_mean_identity(value)
