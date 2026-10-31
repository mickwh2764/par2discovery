"""Tests for the par2 core module."""

import math
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from par2.core import fit_ar2, fit_ar2_batch, classify_dynamics
import pytest

from par2.metrics import half_life, eigenperiod, timescale


def test_unstable_radius_is_not_clipped():
    fit = fit_ar2([3.0 ** i for i in range(12)])
    assert fit["phi1"] == 4
    assert fit["phi2"] == -3
    assert fit["eigenvalue"] == 3
    assert fit["half_life"] is None


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


def test_bootstrap_interval_is_ordered_and_reproducible():
    """Residual bootstrap returns an ordered, seed-reproducible interval."""
    t = np.arange(24)
    x = 0.85 ** t * np.cos(2 * np.pi * t / 12) + 0.1 * np.random.RandomState(3).randn(24)
    a = fit_ar2(x, n_bootstrap=300, seed=1)
    b = fit_ar2(x, n_bootstrap=300, seed=1)
    assert a["eigenvalue_ci"] == b["eigenvalue_ci"]
    assert a["eigenvalue_ci"][0] < a["eigenvalue_ci"][1]
    assert a["phi1_ci"][0] < a["phi1_ci"][1]
    assert a["phi2_ci"][0] < a["phi2_ci"][1]
    assert a["n_timepoints"] == 24


def test_bootstrap_absent_by_default():
    """The interval is opt-in so existing callers are unaffected."""
    x = np.random.RandomState(4).randn(30)
    assert "eigenvalue_ci" not in fit_ar2(x)


def test_bootstrap_interval_narrows_with_series_length():
    """Precision improves with more timepoints; 24 points cannot resolve 0.03."""
    rng = np.random.RandomState(5)

    def series(n):
        t = np.arange(n)
        return np.cos(2 * np.pi * t / 12) * 0.7 ** (t / 12) + 0.3 * rng.randn(n)

    short = fit_ar2(series(24), n_bootstrap=400, seed=2)
    long = fit_ar2(series(96), n_bootstrap=400, seed=2)
    width = lambda r: r["eigenvalue_ci"][1] - r["eigenvalue_ci"][0]
    assert width(long) < width(short)
    assert width(short) > 0.1


def test_timescale_metric():
    assert timescale(math.exp(-1), sampling_interval=2.0) == pytest.approx(2.0)
    assert timescale(0.5, 1.0) * math.log(2) == pytest.approx(half_life(0.5, 1.0))
    assert timescale(0.0) == 0.0
    assert timescale(1.0) == math.inf


def test_timescale_absent_without_sampling_interval():
    x = np.random.RandomState(6).randn(30)
    assert "timescale" not in fit_ar2(x)


def test_timescale_in_fit_matches_metric_and_scales_with_interval():
    x = np.random.RandomState(7).randn(48)
    one = fit_ar2(x, sampling_interval=1.0)
    two = fit_ar2(x, sampling_interval=2.0)
    assert one["timescale"] == pytest.approx(timescale(one["eigenvalue"], 1.0), abs=1e-4)
    assert two["timescale"] == pytest.approx(2 * one["timescale"], abs=1e-3)


def test_timescale_none_for_unstable_fit():
    fit = fit_ar2([3.0 ** i for i in range(12)], sampling_interval=2.0)
    assert fit["timescale"] is None


def test_timescale_interval_and_precise_flag():
    rng = np.random.RandomState(8)
    x = np.zeros(400)
    for t in range(2, 400):
        x[t] = 0.5 * x[t - 1] - 0.3 * x[t - 2] + rng.randn()
    fit = fit_ar2(x, n_bootstrap=300, seed=1, sampling_interval=1.0)
    low, high = fit["timescale_ci"]
    assert low < fit["timescale"] < high
    assert fit["timescale_imprecise"] is False


def test_wide_timescale_interval_warns():
    t = np.arange(12)
    x = np.cos(2 * np.pi * t / 12) + 0.5 * np.random.RandomState(9).randn(12)
    with pytest.warns(UserWarning, match="too short or too noisy"):
        fit = fit_ar2(x, n_bootstrap=300, seed=1, sampling_interval=2.0)
    assert fit["timescale_imprecise"] is True


def test_sampling_interval_must_be_positive():
    with pytest.raises(ValueError, match="sampling_interval"):
        fit_ar2(np.arange(10.0), sampling_interval=0)


def test_batch_carries_timescale():
    m = np.random.RandomState(10).randn(3, 24)
    rows = fit_ar2_batch(m, sampling_interval=2.0)
    assert all("timescale" in r for r in rows)


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
