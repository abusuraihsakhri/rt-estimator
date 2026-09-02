# Rt Estimator

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

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

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`RtEstimate`** — dedicated module for rt estimate evaluation and state verification.

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --task-id <value> --target <value> --primary <value> --secondary <value>
```

### Parameter Reference
- `--task-id`: Specifies input measurement or parameter value.
- `--target`: Specifies input measurement or parameter value.
- `--primary`: Specifies input measurement or parameter value.
- `--secondary`: Specifies input measurement or parameter value.
- `--critical`: Specifies input measurement or parameter value.
- `--status`: Specifies input measurement or parameter value.
- `--input`: Specifies input measurement or parameter value.
- `--output`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `suite_name` | Parameter / observation metric | Required |
| `system_slug` | Parameter / observation metric | Required |
| `standard_reference` | Parameter / observation metric | Required |
| `test_cases` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t rt-estimator .
docker run -p 8000:8000 rt-estimator
```
