# Rt Estimator

Estimate the time-varying effective reproduction number (Rt) from daily incidence using the renewal-equation method of Cori et al. (2013).

This repository includes the scientific Python CLI, a static browser calculator in the web/ directory, and a separate FastAPI task-auditing demonstration retained for compatibility.

## Browser calculator

The browser app uses local JavaScript and requires no backend, Python runtime, API key, external packages or data upload. From the repository root, run:

    python -m http.server 8080 --directory web

Then open http://localhost:8080/. Load the synthetic example or a CSV containing a cases, incidence or count column and, optionally, a date column. Each row represents one consecutive day, including zero-case days. Select serial interval mean and SD, estimation window, Gamma prior and optional trailing moving-average smoothing. Select Calculate Rt to view estimates and download the results.

A GitHub Pages deployment workflow is included in .github/workflows/pages.yml. Configure repository Pages to use GitHub Actions as its deployment source. The live link should only be advertised after the actual deployment has been verified.

## Python installation and use

Requires Python 3.10 or later.

    python -m pip install -r requirements.txt
    python rt_estimator.py --input sample_data/cases.csv --output rt_results.csv
    python rt_estimator.py --input sample_data/cases.csv --si-mean 4.7 --si-sd 2.9 --window 7 --smooth --plot rt_plot.png

Options: --input (required), --output (default rt_estimates.csv), --window (default 7), --si-mean (4.7), --si-sd (2.9), --si-distr (custom discrete serial interval CSV), --prior-mean (5), --prior-sd (5), --smooth, --smooth-window (7), and --plot.

The input CSV must have a header with cases, incidence or count. Counts must be finite and non-negative. Optional dates are copied to the output without checking daily continuity. Custom serial interval weights are supplied one per line beginning with w(0)=0, with positive mass at later lags.

Given daily incidence I(t), the infectiousness Lambda(t) is the sum over positive lags s of w(s) * I(t-s). The Gamma posterior has shape a + sum(I) and rate b + sum(Lambda), with a=(prior mean/prior SD)^2 and b=prior mean/(prior SD)^2. Results include posterior mean and equal-tailed 95% credible intervals.

### Methodological limitations

Estimates assume the specified, fixed serial interval, a renewal model, a Gamma prior and Poisson observations conditional on Rt. Reporting delays, missing observations, imported cases and changes in testing can bias estimates. Smoothing fractionalizes counts and makes the posterior an approximation to the count model. These results are research/educational estimates, not a live surveillance feed or clinical advice.

## Separate audit demonstration

The agents/ and enrichment.py modules implement a generic task-auditing demonstration, not epidemiological Rt estimation. Existing interfaces are retained:

    python cli.py audit --task-id TASK-001 --primary 28.5 --secondary 14.2
    python cli.py batch -i input.csv -o results.csv
    python cli.py verify-audit
    python cli.py serve --host 127.0.0.1 --port 8000

The FastAPI service exposes GET /health, GET /metrics (JSON), POST /api/audit, POST /api/chat, and GET /api/audit/logs.

**Security warning:** The audit API is unauthenticated and not suitable for an untrusted network. Its regex identifier checks do not establish safe de-identification or HIPAA compliance. HMAC-signed audit logs are held in memory; verification detects mutation of retained blocks but cannot prove the completeness of persisted audit history. Do not send patient-identifiable data.

The Docker Compose service binds port 8000 to localhost only and requires AUDIT_SECRET_KEY to be set, for example using a random 32-byte hexadecimal secret, before starting.

## Tests and development

    python -m pytest -q
    node --test web/*.test.mjs
    python -m compileall -q agents rt_estimator.py cli.py enrichment.py simulator.py

CI runs Python 3.10, 3.11 and 3.12. Browser numerical tests use SciPy-computed reference values. The scientific Python CLI requires NumPy and SciPy; Matplotlib is used for optional plotting. The audit demonstration uses FastAPI, Pydantic and Uvicorn.

## Privacy and compatibility

Browser calculations are performed locally: no upload, cookies, analytics, external APIs or persistent browser storage. CSV is downloaded only on request. The app uses modern JavaScript modules supported by current Chrome, Firefox, Safari and Edge. GitHub Pages can host the browser app; the Python CLI and FastAPI server require a Python environment. Pyodide is unnecessary because the browser estimator is implemented in JavaScript.

## Reference

Cori A, Ferguson NM, Fraser C, Cauchemez S. A new framework and software to estimate time-varying reproduction numbers during epidemics. American Journal of Epidemiology. 2013;178(9):1505–1512. doi:10.1093/aje/kwt133.

## License

MIT. See LICENSE.
