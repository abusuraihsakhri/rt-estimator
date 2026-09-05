# Rt Estimator

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## Overview

Estimates the real-time effective reproduction number (Rt) from a daily case
incidence time series, using the renewal-equation method of Cori et al.
(2013), "A New Framework and Software to Estimate Time-Varying Reproduction
Numbers During Epidemics" (American Journal of Epidemiology).

### Method Summary

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

---

## Installation

```bash
# Clone the repository
git clone https://github.com/abusuraihsakhri/rt-estimator.git
cd rt-estimator

# Create and activate a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

---

## Usage

### Core Rt Estimation (CLI)

Estimate Rt from a CSV file containing daily case counts:

```bash
# Basic usage with default parameters
python rt_estimator.py --input sample_data/cases.csv

# With custom serial interval and output file
python rt_estimator.py --input sample_data/cases.csv \
    --output rt_results.csv \
    --si-mean 4.7 --si-sd 2.9 \
    --window 7 \
    --smooth \
    --plot rt_plot.png
```

### CLI Parameters

| Parameter | Description | Default |
|:----------|:------------|:--------|
| `--input` | CSV file with a 'cases' column and optional 'date' column | Required |
| `--output` | Output CSV path | `rt_estimates.csv` |
| `--window` | Sliding window width in days | 7 |
| `--si-mean` | Serial interval mean in days | 4.7 |
| `--si-sd` | Serial interval standard deviation in days | 2.9 |
| `--si-distr` | Path to custom discrete serial interval distribution file | None |
| `--prior-mean` | Gamma prior mean for Rt | 5.0 |
| `--prior-sd` | Gamma prior standard deviation for Rt | 5.0 |
| `--smooth` | Apply 7-day moving average to case counts | False |
| `--smooth-window` | Moving-average window in days | 7 |
| `--plot` | Path to save a PNG plot of incidence and Rt | None |

### Input CSV Format

The input CSV file should contain a header row with at least one of the
following case-count columns (case-insensitive): `cases`, `incidence`, or
`count`. An optional `date` column can be included.

Example:
```csv
date,cases
2026-06-01,5
2026-06-02,6
2026-06-03,8
```

### Enterprise CLI (Multi-Agent System)

The project also includes an enterprise multi-agent system with PHI protection:

```bash
# Run a single task evaluation
python cli.py audit --task-id TASK-001 --primary 28.5 --secondary 14.2

# Batch process CSV records
python cli.py batch -i input.csv -o results.csv

# Verify HMAC audit trail integrity
python cli.py verify-audit

# Launch FastAPI REST server
python cli.py serve --host 127.0.0.1 --port 8000
```

---

## API Server

The FastAPI server exposes the following endpoints:

| Endpoint | Method | Description |
|:---------|:-------|:------------|
| `/health` | GET | Health check |
| `/metrics` | GET | Prometheus-compatible metrics |
| `/api/audit` | POST | Submit a task for evaluation |
| `/api/chat` | POST | Query the supervisory chat |
| `/api/audit/logs` | GET | Retrieve audit trail |

---

## Security Features

* **Zero-PHI Outbound Guard:** Regex-based detection and blocking of SSNs, MRNs, phone numbers, emails, and patient identifiers in outbound data.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Configurable Audit Secret:** Set `AUDIT_SECRET_KEY` environment variable for production deployments.

---

## Testing

Run the automated test suite:

```bash
pytest -v
```

Run the core algorithm tests directly:

```bash
python test_rt_estimator.py
```

---

## Docker Deployment

```bash
# Build the container
docker build -t rt-estimator .

# Run with environment variables
docker run -p 8000:8000 -e AUDIT_SECRET_KEY=your-secret-key rt-estimator

# Or use docker-compose
AUDIT_SECRET_KEY=your-secret-key docker-compose up
```

---

## Project Structure

```
rt-estimator/
├── rt_estimator.py      # Core Rt estimation algorithm
├── cli.py               # Enterprise CLI with multi-agent system
├── test_rt_estimator.py # Core algorithm tests
├── simulator.py         # High-throughput simulation
├── enrichment.py        # Feature enrichment engines
├── requirements.txt     # Python dependencies
├── Dockerfile           # Container build
├── docker-compose.yml   # Container orchestration
├── sample_data/         # Sample case data
├── agents/              # Multi-agent system modules
│   ├── base.py          # Security, PHI guard, audit trail
│   ├── models.py        # Pydantic data models
│   ├── supervisor.py    # Supervisor orchestrator
│   ├── workers.py       # Specialized worker agents
│   ├── api.py           # FastAPI endpoints
│   ├── metrics.py       # Prometheus metrics
│   ├── learning.py      # Bayesian calibration engine
│   ├── llm_factory.py   # LLM provider factory
│   └── streamer.py      # WebSocket telemetry
├── tests/               # Test suite
└── web/                 # Web dashboard
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
