# BLAST Hit Parser & Taxonomic LCA Engine

A pure-Python domain engine and command-line platform for parsing, filtering, translating, and taxonomically classifying NCBI BLAST alignment outputs (Tabular `outfmt 6` / `outfmt 7` / `m8` and XML `outfmt 5`).

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![Tests](https://img.shields.io/badge/Tests-31%20Passing-brightgreen.svg)
![Zero Dependencies](https://img.shields.io/badge/Dependencies-Zero%20External-success.svg)

</div>

---

## 📖 Overview

BLAST (Basic Local Alignment Search Tool) hit parsing is central to pathogen identification, metagenomics, antibiotic resistance profiling, and phylogenetic reconstruction. This repository provides zero-dependency parsing, statistical filtering, sequence translation, and lineage assignment:

- **BLAST Tabular Parsing (`-outfmt 6 / 7` and `m8`):** Parses standard 12-column outputs and custom format strings (TSV and CSV), clustering multiple High-Scoring Segment Pairs (HSPs) per query-subject hit.
- **BLAST XML Parsing (`-outfmt 5`):** Streaming regex-based XML extraction without heavy XML dependencies.
- **Statistical Significance & Bit-Score Transformations:** Full implementation of Karlin-Altschul alignment statistics.
- **HSP Tiling & Overlap Resolution:** Calculates true non-overlapping query coverage across fragmented alignments.
- **Bitscore-Weighted Taxonomic LCA Assignment:** Lowest Common Ancestor resolution with configurable support thresholds and Newick tree export.
- **Six-Frame Translation & ORF Finding:** Translates double-stranded DNA across all reading frames (+1, +2, +3, -1, -2, -3) and extracts open reading frames with start/stop boundary coordinates.

---

## 📐 Bioinformatics & Alignment Statistics

### 1. Karlin-Altschul Expectation Value (E-Value)

The statistical significance of an alignment score $S$ between a query sequence of length $m$ and a database of total length $n$ is evaluated under the Karlin-Altschul extreme value distribution:

$$E = K \cdot m \cdot n \cdot e^{-\lambda S}$$

where:
- $m$: Query sequence length (or effective query length).
- $n$: Total database size in letters (or effective search space size).
- $K$: Scale factor accounting for search space geometry and matrix background frequencies.
- $\lambda$: Natural scale parameter for the scoring system (e.g., BLOSUM62 or nucleotide match/mismatch scores).

### 2. Bit Score Transformation ($S'$)

Raw alignment scores $S$ depend on the scoring matrix and gap penalties. The bit score $S'$ normalizes the raw score into standard information-theoretic units (bits):

$$S' = \frac{\lambda S - \ln K}{\ln 2}$$

Using the bit score, the E-value can be computed independently of scoring matrix parameters:

$$E = m \cdot n \cdot 2^{-S'}$$

### 3. HSP Tiling & Non-Overlapping Query Coverage

When multiple HSPs match a single subject, naive summing of alignment lengths overestimates coverage due to overlapping segment pairs. The parser sorts HSP query intervals $[q_{\text{start}}, q_{\text{end}}]$ and merges overlapping or adjacent segments:

$$\text{Query Coverage} = \frac{\sum_{i=1}^{k} (\text{end}_i - \text{start}_i + 1)}{L_{\text{query}}}$$

where $[\text{start}_i, \text{end}_i]$ are the disjoint merged intervals and $L_{\text{query}}$ is the total query sequence length.

---

## 📋 NCBI BLAST+ Tabular `outfmt 6` Field Schema

Standard NCBI BLAST+ tabular format contains 12 tab- or comma-delimited columns:

| Column | Field Name | Description | Example |
|:------:|:-----------|:-----------------------------------------------------|:-------------------|
| 1 | `qseqid` | Query sequence identifier | `seq1_16S` |
| 2 | `sseqid` | Subject sequence identifier (accession) | `NR_024570.1_Ecoli`|
| 3 | `pident` | Percentage of identical matches (%) | `99.50` |
| 4 | `length` | Alignment length (total base pairs / amino acids) | `1500` |
| 5 | `mismatch` | Number of mismatch positions | `7` |
| 6 | `gapopen` | Number of gap openings | `1` |
| 7 | `qstart` | Start of alignment in query sequence (1-based) | `1` |
| 8 | `qend` | End of alignment in query sequence (1-based) | `1500` |
| 9 | `sstart` | Start of alignment in subject sequence (1-based) | `1` |
| 10 | `send` | End of alignment in subject sequence (1-based) | `1500` |
| 11 | `evalue` | Expectation value (statistical significance) | `0.0` |
| 12 | `bitscore` | Bit score (information content) | `2750.0` |

---

## 💻 CLI Quickstart & Usage

The command-line interface (`cli.py`) provides subcommands for batch processing, tabular parsing, XML extraction, threshold filtering, taxonomic LCA, and six-frame translation.

### 1. Batch Processing & Filtering

Process an input CSV/TSV file of BLAST hits with clinical/bioinformatic quality thresholds:

```bash
python cli.py batch -i sample.csv -o out_filtered.csv --max-evalue 1e-5 --min-bit-score 50.0 --min-identity 70.0
```

### 2. Parse Tabular Alignments (`outfmt 6`)

Inspect hits in standard tabular format or output JSON:

```bash
# Formatted audit table
python cli.py parse-tabular --input sample.csv

# Structured JSON export
python cli.py parse-tabular --input sample.csv --json
```

### 3. Multi-Criteria Hit Filtering

Filter hits by maximum E-value, minimum bit score, percentage identity, and query coverage:

```bash
python cli.py filter --input sample.csv --max-evalue 1e-10 --min-identity 90.0 --min-coverage 0.80
```

### 4. Lowest Common Ancestor (LCA) Taxonomic Assignment

Group hits by query and compute bitscore-weighted taxonomic LCA with Newick phylogenetic tree serialization:

```bash
python cli.py lca --text "seq1\tsubj1\t99.0\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200.0\t100\t100\t\tBacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia"
```

### 5. Six-Frame Translation & ORF Finding

Translate DNA sequences across all 6 frames (+1, +2, +3, -1, -2, -3) and identify open reading frames:

```bash
python cli.py translate --sequence "ATGCGATCGATCGATCGATAGCTAGCTAGCTAATCG" --min-orf-len 5
```

---

## 🐍 Python API Quickstart

```python
from blast_hit_parser import (
    BlastTabularParser,
    FilterCriteria,
    HitFilter,
    KarlinAltschulStatistics,
    SequenceTranslator,
    TaxonomicLCAEngine,
)

# 1. Parse BLAST Tabular Hits
raw_tsv = """
seq1	NR_024570.1_Ecoli	99.50	1500	7	1	1	1500	1	1500	0.0	2750.0
seq1	NR_042817.1_Shigella	99.20	1500	12	0	1	1500	1	1500	0.0	2720.0
seq2	NR_112001.1_Staph	98.00	500	10	0	1	500	1	500	1e-100	920.0
"""
hits = BlastTabularParser.parse_text(raw_tsv)

# 2. Filter Hits by Quality Thresholds
criteria = FilterCriteria(
    max_evalue=1e-10,
    min_bit_score=100.0,
    min_identity=95.0,
    min_query_coverage=0.70,
)
result = HitFilter.filter_hits(hits, criteria)
print(f"Retained {len(result['retained_hits'])} hits ({result['retention_rate_pct']}%)")

# 3. Karlin-Altschul E-value Calculation
evalue = KarlinAltschulStatistics.calculate_evalue(
    bit_score=2750.0,
    query_len=1500,
    db_total_letters=10_000_000,
)
print(f"Calculated E-value: {evalue}")

# 4. Six-Frame Translation
six_frames = SequenceTranslator.translate_six_frames("ATGCCCAAACTGAATTAA")
print("Frame +1:", six_frames[1])
```

---

## 🧪 Testing

Run unit tests across all parsers, filters, statistics, and CLI commands:

```bash
python -m pytest -v -p no:zarr
```

All 31 unit tests run with zero external runtime dependencies.

