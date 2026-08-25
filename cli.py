#!/usr/bin/env python3
"""
Command-Line Interface for BLAST Hit Parser & Taxonomic LCA Engine
===================================================================
Provides comprehensive CLI commands for parsing, filtering, taxonomic LCA
classification, 6-frame translation, and Newick phylogenetic export.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

from blast_hit_parser import (
    BlastTabularParser,
    BlastXMLParser,
    FilterCriteria,
    Hit,
    HitFilter,
    KarlinAltschulStatistics,
    SequenceTranslator,
    TaxonomicLCAEngine,
)


def cmd_parse_tabular(args: argparse.Namespace) -> int:
    if args.input:
        text = Path(args.input).read_text(encoding="utf-8")
    elif args.text:
        text = args.text
    else:
        print("Error: Either --input <file> or --text <raw_blast_text> is required", file=sys.stderr)
        return 1

    hits = BlastTabularParser.parse_text(text)
    if args.json:
        data = [
            {
                "query_id": h.query_id,
                "subject_id": h.subject_id,
                "subject_title": h.subject_title,
                "num_hsps": len(h.hsps),
                "best_evalue": h.best_evalue,
                "total_bit_score": h.total_bit_score,
                "query_coverage": round(h.query_coverage * 100.0, 2),
                "taxonomic_lineage": h.taxonomic_lineage,
            }
            for h in hits
        ]
        print(json.dumps(data, indent=2))
    else:
        print("=" * 75)
        print(f"  BLAST TABULAR PARSER REPORT ({len(hits)} Hits Found)")
        print("=" * 75)
        print(f"{'Query ID':<15} {'Subject ID':<15} {'HSPs':<6} {'E-Value':<12} {'BitScore':<10} {'Cov%':<8}")
        print("-" * 75)
        for h in hits:
            best = h.best_hsp
            e_str = f"{best.evalue:.2e}" if best else "N/A"
            bit_str = f"{best.bit_score:.1f}" if best else "0.0"
            cov_str = f"{h.query_coverage * 100.0:.1f}%"
            print(f"{h.query_id:<15} {h.subject_id:<15} {len(h.hsps):<6} {e_str:<12} {bit_str:<10} {cov_str:<8}")
        print("=" * 75)
    return 0


def cmd_parse_xml(args: argparse.Namespace) -> int:
    if args.input:
        text = Path(args.input).read_text(encoding="utf-8")
    elif args.text:
        text = args.text
    else:
        print("Error: Either --input <file> or --text <raw_xml_text> is required", file=sys.stderr)
        return 1

    hits = BlastXMLParser.parse_xml_text(text)
    if args.json:
        data = [
            {
                "query_id": h.query_id,
                "subject_id": h.subject_id,
                "subject_title": h.subject_title,
                "query_length": h.query_length,
                "subject_length": h.subject_length,
                "hsps": [
                    {
                        "score": hsp.score,
                        "bit_score": hsp.bit_score,
                        "evalue": hsp.evalue,
                        "identity_pct": hsp.identity,
                        "alignment_length": hsp.alignment_length,
                        "gaps": hsp.gaps,
                        "q_start": hsp.q_start,
                        "q_end": hsp.q_end,
                        "s_start": hsp.s_start,
                        "s_end": hsp.s_end,
                    }
                    for hsp in h.hsps
                ],
            }
            for h in hits
        ]
        print(json.dumps(data, indent=2))
    else:
        print("=" * 75)
        print(f"  BLAST XML PARSER REPORT ({len(hits)} Hits Found)")
        print("=" * 75)
        for h in hits:
            print(f"Query: {h.query_id} (len={h.query_length}) -> Subject: {h.subject_id} ({h.subject_title})")
            for i, hsp in enumerate(h.hsps, 1):
                print(f"  HSP {i}: E={hsp.evalue:.2e} | Bits={hsp.bit_score:.1f} | Ident={hsp.identity:.1f}% | "
                      f"AlignLen={hsp.alignment_length} | Gaps={hsp.gaps} | Q:[{hsp.q_start}-{hsp.q_end}] S:[{hsp.s_start}-{hsp.s_end}]")
        print("=" * 75)
    return 0


def cmd_filter(args: argparse.Namespace) -> int:
    if args.input:
        text = Path(args.input).read_text(encoding="utf-8")
    elif args.text:
        text = args.text
    else:
        print("Error: Either --input or --text is required", file=sys.stderr)
        return 1

    hits = BlastTabularParser.parse_text(text)
    criteria = FilterCriteria(
        max_evalue=args.max_evalue,
        min_bit_score=args.min_bit_score,
        min_identity=args.min_identity,
        min_query_coverage=args.min_coverage,
        min_alignment_length=args.min_length,
        max_gaps=args.max_gaps,
    )
    result = HitFilter.filter_hits(hits, criteria)

    if args.json:
        out = {
            "retention_rate_pct": result["retention_rate_pct"],
            "filter_statistics": result["filter_statistics"],
            "retained_hits": [
                {
                    "query_id": h.query_id,
                    "subject_id": h.subject_id,
                    "evalue": h.best_evalue,
                    "bit_score": h.total_bit_score,
                    "identity": h.best_hsp.identity if h.best_hsp else 0.0,
                    "query_coverage": round(h.query_coverage * 100.0, 2),
                }
                for h in result["retained_hits"]
            ],
        }
        print(json.dumps(out, indent=2))
    else:
        stats = result["filter_statistics"]
        print("=" * 60)
        print("  BLAST HIT FILTER AUDIT REPORT")
        print("=" * 60)
        print(f"Total Input Hits        : {stats['total_input_hits']}")
        print(f"Retained Hits           : {stats['retained_hits']} ({result['retention_rate_pct']}%)")
        print("-" * 60)
        print(f"Rejected (E-Value)      : {stats['rejected_evalue']}")
        print(f"Rejected (Bit Score)    : {stats['rejected_bit_score']}")
        print(f"Rejected (% Identity)   : {stats['rejected_identity']}")
        print(f"Rejected (Coverage)     : {stats['rejected_query_coverage']}")
        print(f"Rejected (Min Length)   : {stats['rejected_alignment_length']}")
        print("=" * 60)
    return 0


def cmd_lca(args: argparse.Namespace) -> int:
    if args.input:
        text = Path(args.input).read_text(encoding="utf-8")
    elif args.text:
        text = args.text
    else:
        print("Error: Either --input or --text is required", file=sys.stderr)
        return 1

    hits = BlastTabularParser.parse_text(text)
    # Group by query_id
    by_query: Dict[str, List[Hit]] = {}
    for h in hits:
        by_query.setdefault(h.query_id, []).append(h)

    query_lca_map: Dict[str, Dict[str, Any]] = {}
    for q_id, q_hits in by_query.items():
        query_lca_map[q_id] = TaxonomicLCAEngine.compute_lca(q_hits, min_support_bitscore_fraction=args.min_support)

    newick_str = TaxonomicLCAEngine.export_newick_tree(query_lca_map)

    if args.json:
        print(json.dumps({"assignments": query_lca_map, "newick_tree": newick_str}, indent=2))
    else:
        print("=" * 70)
        print("  TAXONOMIC LOWEST COMMON ANCESTOR (LCA) ASSIGNMENTS")
        print("=" * 70)
        for q_id, res in query_lca_map.items():
            lineage_str = " > ".join(res["lineage_to_lca"]) if res["lineage_to_lca"] else "Unassigned"
            print(f"Query: {q_id}")
            print(f"  LCA Taxon   : {res['lca_taxon']} (Level {res['lca_rank_level']})")
            print(f"  Support     : {res['support_fraction'] * 100:.1f}% (Hits evaluated: {res['hits_evaluated']})")
            print(f"  Lineage Path: {lineage_str}")
            print("-" * 70)
        print("Newick Tree Representation:")
        print(newick_str)
        print("=" * 70)
    return 0


def cmd_translate(args: argparse.Namespace) -> int:
    seq = args.sequence.strip()
    if not seq:
        print("Error: --sequence is required", file=sys.stderr)
        return 1

    six_frames = SequenceTranslator.translate_six_frames(seq)
    orfs = SequenceTranslator.find_orfs(seq, min_length_aa=args.min_orf_len)

    if args.json:
        print(json.dumps({"six_frames": six_frames, "orfs": orfs}, indent=2))
    else:
        print("=" * 65)
        print(f"  SIX-FRAME TRANSLATION & ORF REPORT ({len(seq)} bp)")
        print("=" * 65)
        for f, prot in six_frames.items():
            print(f"Frame {f} ({'+' if f<=3 else '-'}): {prot[:60]}... (len={len(prot)} aa)")
        print("-" * 65)
        print(f"Identified ORFs (min_length >= {args.min_orf_len} aa): {len(orfs)}")
        for i, orf in enumerate(orfs[:5], 1):
            print(f"  ORF #{i} [Frame {orf['frame']}]: {orf['length_aa']} aa ({orf['start_aa']}-{orf['end_aa']}) -> {orf['sequence'][:40]}...")
        print("=" * 65)
    return 0


def cmd_interactive() -> int:
    print("BLAST Hit Parser Interactive CLI")
    print("Commands: tabular, xml, filter, lca, translate, help, exit\n")
    while True:
        try:
            line = input("blast> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nExiting.")
            break
        if not line:
            continue
        if line.lower() in ("exit", "quit"):
            break
        if line.lower() == "help":
            print("Commands: tabular <tsv_line>, translate <dna>, exit")
            continue

        parts = line.split(maxsplit=1)
        cmd = parts[0].lower()
        if cmd == "translate" and len(parts) > 1:
            res = SequenceTranslator.translate_frame(parts[1], 1)
            print(f"Frame 1: {res}")
        elif cmd == "tabular" and len(parts) > 1:
            hits = BlastTabularParser.parse_text(parts[1])
            print(f"Parsed {len(hits)} hits.")
        else:
            print(f"Unrecognized command: {cmd}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="blast_hit_parser_cli",
        description="BLAST Hit Parsing, Filtering & LCA Taxonomic Assignment Platform",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: parse-tabular
    p_tab = subparsers.add_parser("parse-tabular", help="Parse BLAST outfmt 6 / m8 tabular output")
    p_tab.add_argument("--input", "-i", type=str, help="Input tabular file path")
    p_tab.add_argument("--text", "-t", type=str, help="Raw tabular string")
    p_tab.add_argument("--json", action="store_true", help="Output JSON")

    # Subcommand: parse-xml
    p_xml = subparsers.add_parser("parse-xml", help="Parse BLAST outfmt 5 XML output")
    p_xml.add_argument("--input", "-i", type=str, help="Input XML file path")
    p_xml.add_argument("--text", "-t", type=str, help="Raw XML string")
    p_xml.add_argument("--json", action="store_true", help="Output JSON")

    # Subcommand: filter
    p_filt = subparsers.add_parser("filter", help="Filter BLAST hits with configurable thresholds")
    p_filt.add_argument("--input", "-i", type=str, help="Input tabular file path")
    p_filt.add_argument("--text", "-t", type=str, help="Raw tabular string")
    p_filt.add_argument("--max-evalue", type=float, default=1e-5, help="Maximum E-value")
    p_filt.add_argument("--min-bit-score", type=float, default=50.0, help="Minimum bit score")
    p_filt.add_argument("--min-identity", type=float, default=70.0, help="Minimum % identity")
    p_filt.add_argument("--min-coverage", type=float, default=0.50, help="Minimum query coverage")
    p_filt.add_argument("--min-length", type=int, default=30, help="Minimum alignment length")
    p_filt.add_argument("--max-gaps", type=int, default=None, help="Maximum allowed gaps")
    p_filt.add_argument("--json", action="store_true", help="Output JSON")

    # Subcommand: lca
    p_lca = subparsers.add_parser("lca", help="Perform bitscore-weighted LCA taxonomic assignment")
    p_lca.add_argument("--input", "-i", type=str, help="Input tabular file with taxonomy")
    p_lca.add_argument("--text", "-t", type=str, help="Raw tabular text with lineage column")
    p_lca.add_argument("--min-support", type=float, default=0.51, help="Minimum vote fraction for LCA rank")
    p_lca.add_argument("--json", action="store_true", help="Output JSON")

    # Subcommand: translate
    p_trans = subparsers.add_parser("translate", help="Six-frame DNA sequence translation & ORF finder")
    p_trans.add_argument("--sequence", "-s", type=str, required=True, help="DNA nucleotide sequence")
    p_trans.add_argument("--min-orf-len", type=int, default=20, help="Minimum ORF length in amino acids")
    p_trans.add_argument("--json", action="store_true", help="Output JSON")

    # Subcommand: interactive
    subparsers.add_parser("interactive", help="Interactive REPL session")

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "parse-tabular":
        return cmd_parse_tabular(args)
    elif args.command == "parse-xml":
        return cmd_parse_xml(args)
    elif args.command == "filter":
        return cmd_filter(args)
    elif args.command == "lca":
        return cmd_lca(args)
    elif args.command == "translate":
        return cmd_translate(args)
    elif args.command == "interactive":
        return cmd_interactive()
    return 0


if __name__ == "__main__":
    sys.exit(main())
