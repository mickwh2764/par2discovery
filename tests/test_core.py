"""Tests for the par2 core module."""

import math
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from par2.core import fit_ar2, fit_ar2_batch, classify_dynamics
from par2.metrics import half_life, eigenperiod


def test_pure_cosine():
    """A pure cosine with period 12 samples should give |lambda| near 1 and complex roots."""
    t = np.arange(24)
    x = np.cos(2 * np.pi * t / 12)
    result = fit_ar2(x)
    assert result["root_type"] == "Complex"
    assert result["eigenvalue"] > 0.9, f"Expected high eigenvalue, got {result['eigenvalue']}"


def test_white_noise():
    """White noise should give |lambda| near 0."""
    rng = np.random.RandomState(42)
    x = rng.randn(100)
    result = fit_ar2(x)
    assert result["eigenvalue"] < 0.3, f"Expected low eigenvalue for noise, got {result['eigenvalue']}"


def test_damped_oscillation():
    """A damped cosine should give intermediate |lambda| with complex roots."""
    t = np.arange(24)
    x = np.exp(-0.05 * t) * np.cos(2 * np.pi * t / 12) + 0.1 * np.random.RandomState(1).randn(24)
    result = fit_ar2(x)
    assert result["root_type"] == "Complex"
    assert 0.3 < result["eigenvalue"] < 1.0


def test_mean_centering():
    """Result should be identical regardless of a constant offset."""
    t = np.arange(24)
    x = np.cos(2 * np.pi * t / 12)
    r1 = fit_ar2(x)
    r2 = fit_ar2(x + 1000)
    assert abs(r1["eigenvalue"] - r2["eigenvalue"]) < 1e-10


def test_minimum_length():
    """Should raise for fewer than 6 timepoints."""
    try:
        fit_ar2([1, 2, 3, 4, 5])
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


def test_batch():
    """fit_ar2_batch should process multiple genes and sort by eigenvalue."""
    rng = np.random.RandomState(42)
    t = np.arange(24)
    row1 = np.cos(2 * np.pi * t / 12)
    row2 = rng.randn(24)
    matrix = np.vstack([row1, row2])
    results = fit_ar2_batch(matrix, ["oscillator", "noise"])
    assert len(results) == 2
    assert results[0]["eigenvalue"] >= results[1]["eigenvalue"]
    assert results[0]["gene"] == "oscillator"


def test_classify_dynamics():
    assert classify_dynamics(0.95, "Complex") == "Sustained oscillator"
    assert classify_dynamics(0.6, "Complex") == "Damped oscillator"
    assert classify_dynamics(0.6, "Real") == "Overdamped decay"
    assert classify_dynamics(0.2, "Real") == "Rapid decay"
    assert classify_dynamics(1.1, "Complex") == "Unstable"


def test_half_life_metric():
    hl = half_life(0.5, sampling_interval=2.0)
    expected = 2.0 * math.log(2) / math.log(2)
    assert abs(hl - expected) < 0.01


def test_eigenperiod_metric():
    ep = eigenperiod(1.0, -0.5, sampling_interval=2.0)
    assert not math.isnan(ep)
    ep_real = eigenperiod(1.0, 0.5, sampling_interval=2.0)
    assert math.isnan(ep_real)


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for test in tests:
        try:
            test()
            print(f"  PASS  {test.__name__}")
            passed += 1
        except Exception as e:
            print(f"  FAIL  {test.__name__}: {e}")
    print(f"\n{passed}/{len(tests)} tests passed")
