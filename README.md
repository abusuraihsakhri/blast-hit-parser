# Blast Hit Parser

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

BLAST Hit Parser & Taxonomic LCA Engine
========================================
Comprehensive domain engine for parsing, filtering, translating, and taxonomically
classifying BLAST alignment outputs (Tabular outfmt 6 / m8 and XML outfmt 5).

Key Capabilities:
- Standard (12-column) & custom extended BLAST tabular format parsing
- BLAST XML (outfmt 5) parsing without external XML libraries
- Karlin-Altschul E-value statistics & bit-score transformations
- HSP tiling, overlap resolution, and query coverage calculation
- Bitscore-weighted Lowest Common Ancestor (LCA) taxonomic assignment & Newick tree export
- Six-frame translation and Open Reading Frame (ORF) detection

BLAST Hit Parser
Parses BLAST tabular (m8) to best hit, coverage, identity and taxonomic LCA.
Stdlib parser / mapper with batch CSV and single lookup.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`HSP`**: High-Scoring Segment Pair (HSP).
- **`Hit`**: A matched subject sequence hit containing one or more HSPs.
- **`FilterCriteria`**: Filter criteria thresholds for BLAST hits.
- **`BlastTabularParser`**: Parser for BLAST outfmt 6 and m8 tabular reports.
- **`BlastXMLParser`**: Pure-Python regex-based parser for BLAST XML reports (outfmt 5).
- **`KarlinAltschulStatistics`**: Calculates Karlin-Altschul E-values and bit scores.

---

## 📐 Mathematical Formulation & Logic

```text
  """Calculate non-overlapping query coverage fraction (0.0 - 1.0)."""
  bitscore = float(row.get("bitscore", 0.0))
  """Calculates Karlin-Altschul E-values and bit scores."""
  return (lam * raw_score - math.log(k)) / math.log(2)
  score = h.best_hsp.bit_score if h.best_hsp else 1.0
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --input <value> --text <value> --sequence <value> --json <value>
```

### Parameter Reference
- `--input`: Specifies input measurement or parameter value.
- `--text`: Specifies input measurement or parameter value.
- `--sequence`: Specifies input measurement or parameter value.
- `--json`: Specifies input measurement or parameter value.
- `--max-evalue`: Specifies input measurement or parameter value.
- `--min-bit-score`: Specifies input measurement or parameter value.
- `--min-identity`: Specifies input measurement or parameter value.
- `--min-coverage`: Specifies input measurement or parameter value.
- `--min-length`: Specifies input measurement or parameter value.
- `--max-gaps`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `query` | Parameter / observation metric | Required |
| `name` | Parameter / observation metric | Required |

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
docker build -t blast-hit-parser .
docker run -p 8000:8000 blast-hit-parser
```
