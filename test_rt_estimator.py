"""Tests for rt_estimator.py, run with: python test_rt_estimator.py"""

import numpy as np

from rt_estimator import (
    discretize_serial_interval,
    compute_infectiousness,
    moving_average,
    estimate_rt,
)


def test_discretize_serial_interval_is_a_valid_distribution():
    w = discretize_serial_interval(mean=4.7, sd=2.9)
    assert w[0] == 0.0, "zero-day serial interval should have no probability mass"
    assert np.all(w >= 0), "probabilities must be non-negative"
    assert abs(w.sum() - 1.0) < 1e-9, f"distribution must sum to 1, got {w.sum()}"
    # mean of the discretized distribution should be close to the requested mean
    days = np.arange(len(w))
    empirical_mean = np.sum(days * w)
    assert abs(empirical_mean - 4.7) < 0.5, f"empirical mean {empirical_mean} too far from 4.7"


def test_compute_infectiousness_matches_manual_calculation():
    si = np.array([0.0, 0.5, 0.5])  # w(0)=0, w(1)=0.5, w(2)=0.5
    incidence = np.array([10.0, 20.0, 30.0, 40.0])
    lam = compute_infectiousness(incidence, si)

    # lam[0]: no past days -> 0
    assert lam[0] == 0.0
    # lam[1] = w(1)*I[0] = 0.5*10 = 5
    assert abs(lam[1] - 5.0) < 1e-9
    # lam[2] = w(1)*I[1] + w(2)*I[0] = 0.5*20 + 0.5*10 = 15
    assert abs(lam[2] - 15.0) < 1e-9
    # lam[3] = w(1)*I[2] + w(2)*I[1] = 0.5*30 + 0.5*20 = 25
    assert abs(lam[3] - 25.0) < 1e-9


def test_moving_average_trailing_window():
    x = np.array([0.0, 10.0, 20.0, 0.0, 10.0, 20.0])
    ma = moving_average(x, window=3)
    # ma[0] = mean([0]) = 0
    assert abs(ma[0] - 0.0) < 1e-9
    # ma[1] = mean([0,10]) = 5
    assert abs(ma[1] - 5.0) < 1e-9
    # ma[2] = mean([0,10,20]) = 10
    assert abs(ma[2] - 10.0) < 1e-9
    # ma[5] = mean([0,10,20]) = 10
    assert abs(ma[5] - 10.0) < 1e-9


def _simulate_deterministic_renewal(seed_incidence, si, true_r, n_days):
    """Deterministically simulate I_t = true_r * Lambda_t forward in time,
    used to check that estimate_rt recovers a known constant Rt."""
    incidence = list(seed_incidence)
    for _ in range(n_days):
        arr = np.array(incidence, dtype=float)
        lam = compute_infectiousness(arr, si)
        # infectiousness of the *next* day depends on all incidence so far;
        # compute it directly rather than re-running the whole series
        t = len(incidence)
        full = np.array(incidence + [0.0])
        lam_t = compute_infectiousness(full, si)[t]
        next_val = true_r * lam_t
        incidence.append(next_val)
    return np.array(incidence)


def test_estimate_rt_recovers_known_growing_r():
    si = discretize_serial_interval(mean=5.0, sd=2.0, max_days=15)
    true_r = 1.5
    seed = [50.0] * 10
    incidence = _simulate_deterministic_renewal(seed, si, true_r, n_days=40)

    estimates = estimate_rt(incidence, si, window=7, prior_mean=5.0, prior_sd=5.0)
    assert len(estimates) > 0

    # once the transient from the seed/prior has passed, the estimate should
    # converge close to the true R used to generate the data
    last = estimates[-1]
    assert abs(last.mean_r - true_r) < 0.1, f"expected ~{true_r}, got {last.mean_r}"
    # the 95% credible interval should actually contain the true value
    assert last.ci_low <= true_r <= last.ci_high


def test_estimate_rt_recovers_known_declining_r():
    si = discretize_serial_interval(mean=5.0, sd=2.0, max_days=15)
    true_r = 0.7
    seed = [200.0] * 10
    incidence = _simulate_deterministic_renewal(seed, si, true_r, n_days=40)

    estimates = estimate_rt(incidence, si, window=7, prior_mean=5.0, prior_sd=5.0)
    assert len(estimates) > 0

    last = estimates[-1]
    assert abs(last.mean_r - true_r) < 0.1, f"expected ~{true_r}, got {last.mean_r}"
    assert last.ci_low <= true_r <= last.ci_high


def test_estimate_rt_stable_incidence_gives_r_near_one():
    # a serial-interval-weighted sum over a constant incidence series should
    # approximately equal the incidence itself (weights sum to 1), so Rt
    # should be close to 1
    si = discretize_serial_interval(mean=4.7, sd=2.9, max_days=15)
    incidence = np.array([100.0] * 30)
    estimates = estimate_rt(incidence, si, window=7, prior_mean=5.0, prior_sd=5.0)
    last = estimates[-1]
    assert abs(last.mean_r - 1.0) < 0.05, f"expected Rt~1.0 for stable incidence, got {last.mean_r}"


def test_estimate_rt_rejects_window_larger_than_series():
    si = np.array([0.0, 1.0])
    try:
        estimate_rt([1.0, 2.0, 3.0], si, window=10)
        assert False, "expected ValueError for window longer than series"
    except ValueError:
        pass


if __name__ == "__main__":
    tests = [obj for name, obj in list(globals().items()) if name.startswith("test_")]
    failures = 0
    for t in tests:
        try:
            t()
            print(f"PASS: {t.__name__}")
        except AssertionError as e:
            failures += 1
            print(f"FAIL: {t.__name__}: {e}")
    if failures:
        raise SystemExit(f"{failures} test(s) failed")
    print(f"All {len(tests)} tests passed.")
