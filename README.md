# BLAST Hit Parser & Taxonomic LCA Engine

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python: 3.10+](https://img.shields.io/badge/Python-3.10%2B-brightgreen.svg)](https://python.org)
[![Tests: 30 Passing](https://img.shields.io/badge/Tests-30%20Passing-success.svg)](test_blast_hit_parser.py)
[![Domain: Bioinformatics & Genomics](https://img.shields.io/badge/Domain-Bioinformatics%20%26%20Metagenomics-blueviolet.svg)](#)

A high-performance parser, statistical evaluator, and taxonomic classifier for NCBI BLAST alignments supporting both **Tabular (outfmt 6 / m8)** and **XML (outfmt 5)** formats.

---

## Key Features

1. **BLAST Tabular Parsing (`-outfmt 6` / `m8`)**
   - Standard 12-column parsing (`qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore`).
   - Extended and custom field parsing (`qlen, slen, stitle, lineage, taxid`).
   - Automated grouping of multiple High-Scoring Pairs (HSPs) per query-subject hit.

2. **BLAST XML Parsing (`-outfmt 5`)**
   - Zero-dependency XML iteration and hit extraction without requiring `lxml`.
   - Comprehensive extraction of HSP coordinates, bit scores, and E-values.

3. **Karlin-Altschul E-Value Statistics & Query Coverage**
   - Mathematical transformation between raw alignment score $S$ and bit score $S'$:
     $$S' = \frac{\lambda S - \ln K}{\ln 2}$$
   - Expectation value computation $E = m \cdot n \cdot 2^{-S'}$.
   - Non-overlapping interval union algorithm for accurate query coverage percentage calculation across multiple HSPs.

4. **Bitscore-Weighted Lowest Common Ancestor (LCA)**
   - Rank-by-rank weighted voting algorithm across taxonomic lineages.
   - Configurable minimum support fraction thresholding for high-confidence metagenomic assignments.
   - Newick phylogenetic tree serialization from assigned lineage paths.

5. **Six-Frame Translation & Open Reading Frame (ORF) Detection**
   - Translation across all 6 reading frames (frames 1–3 forward, frames 4–6 reverse complement).
   - ORF extraction with configurable minimum peptide length and start codon constraints.

---

## Installation

```bash
git clone https://github.com/example/blast-hit-parser.git
cd blast-hit-parser
```

*Requires Python 3.10+ with zero external third-party dependencies (pure standard library).*

---

## Command-Line Interface (CLI)

```bash
# 1. Parse BLAST tabular file
python cli.py parse-tabular --input results.m8

# 2. Parse BLAST XML file and output structured JSON
python cli.py parse-xml --input blast_results.xml --json

# 3. Filter hits by E-value, bit score, identity, and query coverage
python cli.py filter --input results.m8 --max-evalue 1e-10 --min-bit-score 100 --min-identity 90.0

# 4. Perform taxonomic Lowest Common Ancestor (LCA) assignment and export Newick tree
python cli.py lca --input results_with_taxonomy.tsv --min-support 0.60

# 5. Translate nucleotide sequence in 6 frames and identify ORFs
python cli.py translate --sequence "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG" --min-orf-len 5
```

---

## Python API Usage

```python
from blast_hit_parser import (
    BlastTabularParser,
    BlastXMLParser,
    HitFilter,
    FilterCriteria,
    TaxonomicLCAEngine,
    SequenceTranslator,
)

# Parse Tabular BLAST output
tabular_data = """
q1	NR_024570.1_Ecoli	99.50	1500	7	1	1	1500	1	1500	0.0	2750.0
q1	NR_042817.1_Shigella	99.20	1500	12	0	1	1500	1	1500	0.0	2720.0
"""
hits = BlastTabularParser.parse_text(tabular_data)
for hit in hits:
    print(f"Hit: {hit.subject_id} | Best E-value: {hit.best_evalue} | Coverage: {hit.query_coverage*100:.1f}%")

# Filter hits with strict criteria
criteria = FilterCriteria(max_evalue=1e-50, min_identity=95.0, min_bit_score=500.0)
filtered = HitFilter.filter_hits(hits, criteria)
print(f"Retained {len(filtered['retained_hits'])} of {len(hits)} hits")

# Six-frame translation
frames = SequenceTranslator.translate_six_frames("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
print("Frame 1 peptide:", frames[1])
```

---

## Test Suite

Run the full unit test suite:

```bash
python -m unittest test_blast_hit_parser.py
```

All 30 unit tests verify:
- Tabular 12-column and XML outfmt 5 parsing accuracy
- HSP interval merging and non-overlapping query coverage
- Karlin-Altschul mathematical relationships
- Multi-criteria rejection logging
- Bitscore-weighted LCA voting and Newick string construction
- 6-frame translation and ORF detection

---

## License

MIT License. See `LICENSE` for details.
