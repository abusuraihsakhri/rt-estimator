"""
Rt Estimator
============

Estimates the real-time effective reproduction number (Rt) from a daily case
incidence time series, using the renewal-equation method of Cori et al.
(2013), "A New Framework and Software to Estimate Time-Varying Reproduction
Numbers During Epidemics" (American Journal of Epidemiology).

Method summary
--------------
Given daily incidence I_t and a discrete serial interval distribution w(s)
(the probability that the interval between the onset of a primary and a
secondary case is s days), the total infectiousness at time t is

    Lambda_t = sum_{s=1}^{W} w(s) * I_{t-s}

Assuming case counts in a sliding window of tau days follow a Poisson
process with rate R * Lambda_t, and placing a Gamma(a, b) prior on R
(shape a, rate b), the posterior distribution of R over the window
[t - tau + 1, t] is exactly

    R | window ~ Gamma(a + sum(I), 1 / (1/b + sum(Lambda)))

This module implements that estimator, a standard discretization of a
Gamma-distributed serial interval, optional 7-day moving-average smoothing
of noisy case counts, and CSV/plot output.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass

import numpy as np
from scipy.stats import gamma as gamma_dist


# ---------------------------------------------------------------------------
# Serial interval
# ---------------------------------------------------------------------------

def discretize_serial_interval(mean: float, sd: float, max_days: int | None = None) -> np.ndarray:
    """Discretize a Gamma-distributed serial interval.

    Uses the standard discretization (as implemented in the EpiEstim R
    package's DiscrSI function): for a continuous Gamma(shape, scale)
    distribution with the given mean and standard deviation, the probability
    mass assigned to an integer day k is

        w(k) = k*F(k) + (k-2)*F(k-2) - 2*(k-1)*F(k-1)
               + shape*scale*(2*F1(k-1) - F1(k-2) - F1(k))

    where F is the CDF of Gamma(shape, scale), F1 is the CDF of
    Gamma(shape+1, scale), and F(x) = 0 for x < 0.

    Returns an array w of length max_days+1 where w[0] == 0 (no
    zero-day serial interval) and w[1:] sums to 1.
    """
    if mean <= 1:
        raise ValueError("serial interval mean must be > 1 day")
    if sd <= 0:
        raise ValueError("serial interval sd must be > 0")

    shape = ((mean - 1) / sd) ** 2
    scale = (sd ** 2) / (mean - 1)

    if max_days is None:
        max_days = int(np.ceil(mean + 10 * sd))
    max_days = max(int(max_days), 1)

    ks = np.arange(0, max_days + 1, dtype=float)

    def cdf(k, a):
        k = np.asarray(k, dtype=float)
        out = np.zeros_like(k)
        mask = k >= 0
        out[mask] = gamma_dist.cdf(k[mask], a=a, scale=scale)
        return out

    w = ks * cdf(ks, shape) + (ks - 2) * cdf(ks - 2, shape) - 2 * (ks - 1) * cdf(ks - 1, shape)
    w += shape * scale * (2 * cdf(ks - 1, shape + 1) - cdf(ks - 2, shape + 1) - cdf(ks, shape + 1))
    w = np.maximum(w, 0.0)
    w[0] = 0.0

    total = w.sum()
    if total <= 0:
        raise ValueError("serial interval discretization failed for the given mean/sd")
    return w / total


def load_serial_interval_distribution(path: str) -> np.ndarray:
    """Load a custom discrete serial interval distribution from a file.

    The file must contain one non-negative number per line (or per row of a
    single-column CSV), giving w(0), w(1), w(2), ... The array is
    normalized so it sums to 1.
    """
    import os
    if not os.path.isfile(path):
        raise FileNotFoundError(f"serial interval file not found: {path}")
    values = []
    with open(path, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            cell = row[0].strip()
            if cell == "" or not _is_float(cell):
                continue
            values.append(float(cell))
    if not values:
        raise ValueError(f"no numeric values found in serial interval file: {path}")
    w = np.array(values, dtype=float)
    if np.any(w < 0):
        raise ValueError("serial interval distribution values must be non-negative")
    total = w.sum()
    if total <= 0:
        raise ValueError("serial interval distribution must sum to a positive value")
    return w / total


def _is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


# ---------------------------------------------------------------------------
# Smoothing
# ---------------------------------------------------------------------------

def moving_average(x, window: int = 7) -> np.ndarray:
    """Trailing (causal) moving average: out[t] = mean(x[max(0,t-window+1):t+1]).

    Trailing rather than centered so that it can be applied to real-time
    data without looking into the future.
    """
    x = np.asarray(x, dtype=float)
    if window < 1:
        raise ValueError("window must be >= 1")
    out = np.zeros_like(x)
    for t in range(len(x)):
        start = max(0, t - window + 1)
        out[t] = x[start:t + 1].mean()
    return out


# ---------------------------------------------------------------------------
# Renewal equation / Rt estimation
# ---------------------------------------------------------------------------

def compute_infectiousness(incidence, si_distr) -> np.ndarray:
    """Compute Lambda_t = sum_{s=1}^{W} w(s) * I_{t-s} for every t.

    si_distr[0] is assumed to be 0 (or is ignored); si_distr[s] is the
    probability the serial interval equals s days.
    """
    incidence = np.asarray(incidence, dtype=float)
    si = np.asarray(si_distr, dtype=float)
    T = len(incidence)
    W = len(si)
    lam = np.zeros(T)
    for t in range(T):
        s_max = min(t, W - 1)
        if s_max >= 1:
            s = np.arange(1, s_max + 1)
            lam[t] = np.sum(si[s] * incidence[t - s])
    return lam


@dataclass
class RtEstimate:
    t_start: int
    t_end: int
    mean_r: float
    ci_low: float
    ci_high: float
    sum_incidence: float
    sum_lambda: float


def estimate_rt(
    incidence,
    si_distr,
    window: int = 7,
    prior_mean: float = 5.0,
    prior_sd: float = 5.0,
) -> list[RtEstimate]:
    """Estimate Rt over sliding windows using the Cori et al. (2013) method.

    For each window [t_end - window + 1, t_end], the posterior of R is
    Gamma(a + sum(I), 1 / (1/b + sum(Lambda))), where the Gamma(a, b) prior
    (shape a, rate b) is parameterized here by its mean and standard
    deviation. Windows where the sum of infectiousness is zero (e.g. before
    any cases have occurred) are skipped, since R is not identifiable there.
    """
    incidence = np.asarray(incidence, dtype=float)
    T = len(incidence)
    if window < 1:
        raise ValueError("window must be >= 1")
    if T < window:
        raise ValueError(
            f"incidence series has only {T} points, shorter than window={window}"
        )
    if prior_mean <= 0 or prior_sd <= 0:
        raise ValueError("prior mean and sd must be positive")

    lam = compute_infectiousness(incidence, si_distr)

    # Gamma(shape, rate) prior parameterized by mean/sd:
    #   mean = shape / rate, var = shape / rate^2
    shape_prior = (prior_mean / prior_sd) ** 2
    rate_prior = prior_mean / (prior_sd ** 2)

    results: list[RtEstimate] = []
    for t_end in range(window - 1, T):
        t_start = t_end - window + 1
        sum_I = float(incidence[t_start:t_end + 1].sum())
        sum_lambda = float(lam[t_start:t_end + 1].sum())
        if sum_lambda <= 0:
            continue

        post_shape = shape_prior + sum_I
        post_rate = rate_prior + sum_lambda
        scale = 1.0 / post_rate

        mean_r = post_shape * scale
        ci_low = float(gamma_dist.ppf(0.025, post_shape, scale=scale))
        ci_high = float(gamma_dist.ppf(0.975, post_shape, scale=scale))

        results.append(
            RtEstimate(
                t_start=t_start,
                t_end=t_end,
                mean_r=mean_r,
                ci_low=ci_low,
                ci_high=ci_high,
                sum_incidence=sum_I,
                sum_lambda=sum_lambda,
            )
        )
    return results


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_incidence_csv(path: str):
    """Load a daily case incidence series from a CSV file.

    Expects a header row containing a case-count column named (case
    insensitive) one of "cases", "incidence", or "count", and optionally a
    "date" column. If no date column is present, dates are returned as
    None and day indices are used instead.
    """
    import os
    if not os.path.isfile(path):
        raise FileNotFoundError(f"incidence CSV file not found: {path}")
    dates = []
    cases = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"could not read a header row from {path}")
        lower_map = {name.lower().strip(): name for name in reader.fieldnames}

        case_col = None
        for candidate in ("cases", "incidence", "count"):
            if candidate in lower_map:
                case_col = lower_map[candidate]
                break
        if case_col is None:
            raise ValueError(
                "input CSV must have a column named 'cases', 'incidence', or 'count'"
            )
        date_col = lower_map.get("date")

        for row in reader:
            raw = row[case_col]
            if raw is None or raw.strip() == "":
                continue
            cases.append(float(raw))
            dates.append(row[date_col].strip() if date_col else None)

    if not cases:
        raise ValueError(f"no case count rows found in {path}")
    return dates, np.array(cases, dtype=float)


def write_rt_csv(path: str, estimates: list[RtEstimate], dates=None):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["t_start", "t_end", "date", "mean_r", "ci_low_95", "ci_high_95",
             "sum_incidence", "sum_infectiousness"]
        )
        for e in estimates:
            date_str = ""
            if dates is not None and e.t_end < len(dates) and dates[e.t_end] is not None:
                date_str = dates[e.t_end]
            writer.writerow(
                [
                    e.t_start,
                    e.t_end,
                    date_str,
                    f"{e.mean_r:.6f}",
                    f"{e.ci_low:.6f}",
                    f"{e.ci_high:.6f}",
                    e.sum_incidence,
                    f"{e.sum_lambda:.6f}",
                ]
            )


def plot_rt(incidence, estimates: list[RtEstimate], out_path: str,
            smoothed_incidence=None, dates=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    x_incidence = np.arange(len(incidence))
    ax1.bar(x_incidence, incidence, color="#8fb8de", label="daily cases")
    if smoothed_incidence is not None:
        ax1.plot(x_incidence, smoothed_incidence, color="#1f4e79", linewidth=2,
                  label="7-day moving average")
    ax1.set_ylabel("daily case incidence")
    ax1.legend(loc="upper left")

    t_end = np.array([e.t_end for e in estimates])
    mean_r = np.array([e.mean_r for e in estimates])
    ci_low = np.array([e.ci_low for e in estimates])
    ci_high = np.array([e.ci_high for e in estimates])

    ax2.plot(t_end, mean_r, color="#c0392b", linewidth=2, label="Rt (mean)")
    ax2.fill_between(t_end, ci_low, ci_high, color="#c0392b", alpha=0.2,
                      label="95% credible interval")
    ax2.axhline(1.0, color="black", linestyle="--", linewidth=1)
    ax2.set_ylabel("effective reproduction number (Rt)")
    ax2.set_xlabel("day index")
    ax2.legend(loc="upper left")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="rt_estimator",
        description=(
            "Estimate the real-time effective reproduction number (Rt) from a "
            "daily case incidence time series using the Cori et al. (2013) "
            "renewal-equation method."
        ),
    )
    p.add_argument("--input", required=True, help="CSV file with a 'cases' column and optional 'date' column")
    p.add_argument("--output", default="rt_estimates.csv", help="output CSV path (default: rt_estimates.csv)")
    p.add_argument("--window", type=int, default=7, help="sliding window width in days (default: 7)")
    p.add_argument("--si-mean", type=float, default=4.7,
                    help="serial interval mean in days (default: 4.7, a commonly used COVID-19 estimate)")
    p.add_argument("--si-sd", type=float, default=2.9,
                    help="serial interval standard deviation in days (default: 2.9)")
    p.add_argument("--si-distr", default=None,
                    help="path to a file with a custom discrete serial interval distribution "
                         "(one probability per line, w(0), w(1), ...); overrides --si-mean/--si-sd")
    p.add_argument("--prior-mean", type=float, default=5.0, help="Gamma prior mean for Rt (default: 5.0)")
    p.add_argument("--prior-sd", type=float, default=5.0, help="Gamma prior standard deviation for Rt (default: 5.0)")
    p.add_argument("--smooth", action="store_true", help="apply a 7-day moving average to case counts before estimation")
    p.add_argument("--smooth-window", type=int, default=7, help="moving-average window in days (default: 7)")
    p.add_argument("--plot", default=None, help="optional path to save a PNG plot of incidence and Rt")
    return p


def main(argv=None):
    args = build_arg_parser().parse_args(argv)

    dates, incidence = load_incidence_csv(args.input)

    if args.si_distr:
        si = load_serial_interval_distribution(args.si_distr)
    else:
        si = discretize_serial_interval(args.si_mean, args.si_sd)

    smoothed = None
    estimation_input = incidence
    if args.smooth:
        smoothed = moving_average(incidence, window=args.smooth_window)
        estimation_input = smoothed

    try:
        estimates = estimate_rt(
            estimation_input,
            si,
            window=args.window,
            prior_mean=args.prior_mean,
            prior_sd=args.prior_sd,
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if not estimates:
        print("error: no Rt estimates could be produced (check that cases are non-zero)", file=sys.stderr)
        return 1

    write_rt_csv(args.output, estimates, dates=dates)
    print(f"wrote {len(estimates)} Rt estimates to {args.output}")

    last = estimates[-1]
    print(
        f"latest estimate (day {last.t_end}): "
        f"Rt = {last.mean_r:.2f} (95% CI {last.ci_low:.2f}-{last.ci_high:.2f})"
    )

    if args.plot:
        plot_rt(incidence, estimates, args.plot, smoothed_incidence=smoothed, dates=dates)
        print(f"wrote plot to {args.plot}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
