#!/usr/bin/env python3
"""
BLAST Hit Filtering & LCA Taxonomic Assignment
Threshold pipeline (e-value, bit-score, query coverage, identity) with
per-stage removal stats, bitscore-weighted lowest-common-ancestor assignment,
and Newick serialization of the resulting taxonomic tree.

Zero-dependency. Author: Dr. Abu Suraih Sakhri. License: MIT.
"""
import argparse
import json
import math
import sys
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Tuple


@dataclass
class BlastHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    qlen: int
    evalue: float
    bitscore: float
    lineage: List[str] = field(default_factory=list)   # root->...->species

    @property
    def qcov(self) -> float:
        return self.length / self.qlen if self.qlen else 0.0


@dataclass
class FilterConfig:
    min_evalue: float = 1e-5
    min_bitscore: float = 50.0
    min_qcov: float = 0.70
    min_pident: float = 85.0

    def strict(self) -> "FilterConfig":
        return FilterConfig(min_evalue=1e-20, min_bitscore=100.0,
                            min_qcov=0.90, min_pident=95.0)


def filter_hits(hits: List[BlastHit], cfg: FilterConfig = None,
                strict: bool = False) -> Dict[str, Any]:
    cfg = (cfg.strict() if strict and cfg else cfg or FilterConfig())
    if strict:
        cfg = cfg.strict()
    stats = {"input": len(hits), "removed_evalue": 0, "removed_bitscore": 0,
             "removed_coverage": 0, "removed_identity": 0}
    kept = []
    for h in hits:
        if h.evalue > cfg.min_evalue:
            stats["removed_evalue"] += 1; continue
        if h.bitscore < cfg.min_bitscore:
            stats["removed_bitscore"] += 1; continue
        if h.qcov < cfg.min_qcov:
            stats["removed_coverage"] += 1; continue
        if h.pident < cfg.min_pident:
            stats["removed_identity"] += 1; continue
        kept.append(h)
    stats["kept"] = len(kept)
    return {"config": vars(cfg), "stats": stats, "kept": kept}


def lca_assign(hits: List[BlastHit], min_support_bitscore_fraction: float = 0.05
               ) -> Dict[str, Any]:
    """Bitscore-weighted LCA across a query's hits using rank-by-rank voting."""
    if not hits:
        return {"lca": None, "note": "no hits passed filtering"}
    total_bits = sum(h.bitscore for h in hits)
    depth = max(len(h.lineage) for h in hits)
    assigned_ranks = {}
    for lvl in range(depth):
        votes: Dict[str, float] = {}
        for h in hits:
            if lvl >= len(h.lineage):
                break   # shallower lineage cannot vote deeper than itself
            taxon = h.lineage[lvl]
            votes[taxon] = votes.get(taxon, 0.0) + h.bitscore
        best_taxon, best_score = None, 0.0
        for t, s in votes.items():
            if s > best_score and s / total_bits >= min_support_bitscore_fraction:
                best_taxon, best_score = t, s
        if best_taxon is None:
            break
        assigned_ranks[lvl] = best_taxon
    lca_rank_idx = max(assigned_ranks.keys(), default=-1)
    return {
        "lca": assigned_ranks.get(lca_rank_idx),
        "lineage_to_lca": [assigned_ranks[i] for i in sorted(assigned_ranks)],
        "support_fraction": round(
            sum(h.bitscore for h in hits
                if lca_rank_idx < len(h.lineage) and
                h.lineage[lca_rank_idx] == assigned_ranks.get(lca_rank_idx)) / max(total_bits, 1e-9), 4),
        "hits_used": len(hits),
    }


def build_tree_and_newick(assignments: Dict[str, Dict[str, Any]]) -> str:
    """Nested-dict tree from per-query lineages -> Newick string."""
    tree: Dict[str, Any] = {}
    for q, res in assignments.items():
        node = tree
        for taxon in res["lineage_to_lca"]:
            node = node.setdefault(taxon, {})
    def render(node: Dict[str, Any]) -> str:
        if not node:
            return ""
        children = ",".join(f"{c}{render(sub)}" for c, sub in node.items())
        return f"({children})"
    inner = render(tree).rstrip()
    leaves = []

    def collect(node, path):
        if not node:
            leaves.append(path[-1])
            return
        for c, sub in node.items():
            collect(sub, path + [c])
    collect(tree, [])
    return f"{inner};".replace("()", "") if tree else ";"


def best_hit_per_query(hits: List[BlastHit]) -> Dict[str, BlastHit]:
    out: Dict[str, BlastHit] = {}
    for h in hits:
        cur = out.get(h.qseqid)
        if cur is None or h.bitscore > cur.bitscore:
            out[h.qseqid] = h
    return out


def summarize(hits: List[BlastHit]) -> Dict[str, Any]:
    if not hits:
        return {"queries": 0}
    bh = best_hit_per_query(hits)
    return {
        "queries": len(bh),
        "mean_hits_per_query": round(len(hits) / len(bh), 2),
        "mean_best_bitscore": round(sum(h.bitscore for h in bh.values()) / len(bh), 1),
        "mean_best_identity": round(sum(h.pident for h in bh.values()) / len(bh), 2),
    }


if __name__ == "__main__":
    demo_hits = [
        BlastHit("q1", "sp|P1|Ecoli", 99.2, 480, 480, 0.0, 940,
                 ["root", "Bacteria", "Proteobacteria", "Gammaproteobacteria",
                  "Enterobacterales", "Escherichia"]),
        BlastHit("q1", "sp|P2|Shigella", 98.8, 476, 480, 0.0, 930,
                 ["root", "Bacteria", "Proteobacteria", "Gammaproteobacteria",
                  "Enterobacterales", "Shigella"]),
        BlastHit("q1", "sp|P3|Salmonella", 91.4, 460, 480, 1e-160, 700,
                 ["root", "Bacteria", "Proteobacteria", "Gammaproteobacteria",
                  "Enterobacterales", "Salmonella"]),
        BlastHit("q2", "sp|P4|StaphAureus", 97.0, 300, 450, 0.0, 520,
                 ["root", "Bacteria", "Firmicutes", "Bacilli", "Staphylococcus"]),
        BlastHit("q2", "sp|P5|LowCov", 60.0, 120, 450, 0.02, 45,
                 ["root", "Bacteria"]),
    ]

    filtered = filter_hits(demo_hits, strict=True)
    print(json.dumps(filtered["stats"], indent=2))

    assignments = {}
    by_query: Dict[str, List[BlastHit]] = {}
    for h in filtered["kept"]:
        by_query.setdefault(h.qseqid, []).append(h)
    for q, hs in by_query.items():
        assignments[q] = lca_assign(hs)
    print(json.dumps(assignments, indent=2))
    print(json.dumps({"newick": build_tree_and_newick(assignments)}, indent=2))
    print(json.dumps(summarize(filtered["kept"]), indent=2))
