"""
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
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Any, Dict, Iterator, List, Optional, Tuple


# Standard NCBI Genetic Code (Codon Table 1)
STANDARD_CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

DNA_COMPLEMENT: Dict[str, str] = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "a": "t", "t": "a", "g": "c", "c": "g",
    "N": "N", "n": "n", "U": "A", "u": "a",
}


@dataclass
class HSP:
    """High-Scoring Segment Pair (HSP)."""
    score: float
    bit_score: float
    evalue: float
    identity: float  # percentage, e.g., 98.5
    alignment_length: int
    mismatches: int = 0
    gaps: int = 0
    q_start: int = 1
    q_end: int = 1
    s_start: int = 1
    s_end: int = 1
    q_frame: int = 1
    s_frame: int = 1


@dataclass
class Hit:
    """A matched subject sequence hit containing one or more HSPs."""
    query_id: str
    subject_id: str
    subject_title: str = ""
    query_length: int = 0
    subject_length: int = 0
    hsps: List[HSP] = field(default_factory=list)
    taxonomic_lineage: List[str] = field(default_factory=list)

    @property
    def best_hsp(self) -> Optional[HSP]:
        return min(self.hsps, key=lambda h: (h.evalue, -h.bit_score)) if self.hsps else None

    @property
    def best_evalue(self) -> float:
        h = self.best_hsp
        return h.evalue if h else float("inf")

    @property
    def total_bit_score(self) -> float:
        return sum(h.bit_score for h in self.hsps)

    @property
    def query_coverage(self) -> float:
        """Calculate non-overlapping query coverage fraction (0.0 - 1.0)."""
        if self.query_length <= 0 or not self.hsps:
            return 0.0
        intervals = sorted([(min(h.q_start, h.q_end), max(h.q_start, h.q_end)) for h in self.hsps])
        merged = []
        for start, end in intervals:
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        covered_bases = sum(end - start + 1 for start, end in merged)
        return min(1.0, covered_bases / self.query_length)


@dataclass
class FilterCriteria:
    """Filter criteria thresholds for BLAST hits."""
    max_evalue: float = 1e-5
    min_bit_score: float = 50.0
    min_identity: float = 70.0
    min_query_coverage: float = 0.50
    min_alignment_length: int = 30
    max_gaps: Optional[int] = None


# ==============================================================================
# 1. TABULAR PARSER (-outfmt 6 / m8)
# ==============================================================================

class BlastTabularParser:
    """Parser for BLAST outfmt 6 and m8 tabular reports."""

    STANDARD_FIELDS = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]

    @classmethod
    def parse_text(
        cls,
        tabular_text: str,
        custom_fields: Optional[List[str]] = None,
        delimiter: str = "\t",
    ) -> List[Hit]:
        """
        Parse raw tabular BLAST text into a list of Hit objects grouped by query and subject.
        """
        hits_map: Dict[Tuple[str, str], Hit] = {}
        fields = custom_fields or cls.STANDARD_FIELDS

        for line in tabular_text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                # Handle comment line defining custom fields if present
                if line.startswith("# Fields:"):
                    fields = [f.strip() for f in line.replace("# Fields:", "").split(",")]
                continue

            # Support comma-separated (CSV), tab-separated (TSV), or whitespace-separated lines
            if delimiter == "\t" and "\t" not in line and "," in line:
                parts = [p.strip() for p in line.split(",")]
            else:
                parts = [p.strip() for p in line.split(delimiter)]
            if len(parts) < len(fields):
                # Try whitespace delimiter fallback
                parts = line.split()
            if len(parts) < 12 and not custom_fields:
                continue

            row = dict(zip(fields, parts))

            # Skip header line if present (e.g. qseqid, sseqid, etc.)
            if row.get("qseqid", "").lower() in ("qseqid", "query", "query_id") and row.get("pident", "").lower() in ("pident", "% identity", "identity"):
                continue

            qseqid = row.get("qseqid", "unknown_query")
            sseqid = row.get("sseqid", "unknown_subject")
            pident = float(row.get("pident", 0.0))
            length = int(row.get("length", 0))
            mismatch = int(row.get("mismatch", 0))
            gapopen = int(row.get("gapopen", 0))
            qstart = int(row.get("qstart", 1))
            qend = int(row.get("qend", length))
            sstart = int(row.get("sstart", 1))
            send = int(row.get("send", length))
            evalue = float(row.get("evalue", 1.0))
            bitscore = float(row.get("bitscore", 0.0))
            qlen = int(row.get("qlen", max(qstart, qend)))
            slen = int(row.get("slen", max(sstart, send)))
            stitle = row.get("stitle", "")

            # Optional taxonomy fields
            tax_str = row.get("sscinames", row.get("lineage", ""))
            lineage = [t.strip() for t in tax_str.split(";") if t.strip()] if tax_str else []

            hsp = HSP(
                score=bitscore,
                bit_score=bitscore,
                evalue=evalue,
                identity=pident,
                alignment_length=length,
                mismatches=mismatch,
                gaps=gapopen,
                q_start=qstart,
                q_end=qend,
                s_start=sstart,
                s_end=send,
            )

            key = (qseqid, sseqid)
            if key not in hits_map:
                hits_map[key] = Hit(
                    query_id=qseqid,
                    subject_id=sseqid,
                    subject_title=stitle,
                    query_length=qlen,
                    subject_length=slen,
                    hsps=[hsp],
                    taxonomic_lineage=lineage,
                )
            else:
                hit = hits_map[key]
                hit.hsps.append(hsp)
                if qlen > hit.query_length:
                    hit.query_length = qlen
                if slen > hit.subject_length:
                    hit.subject_length = slen

        return list(hits_map.values())


# ==============================================================================
# 2. XML PARSER (-outfmt 5)
# ==============================================================================

class BlastXMLParser:
    """Pure-Python regex-based parser for BLAST XML reports (outfmt 5)."""

    @classmethod
    def parse_xml_text(cls, xml_text: str) -> List[Hit]:
        """Parse XML output into list of Hit objects."""
        hits: List[Hit] = []

        # Find iterations (per query)
        iteration_blocks = re.findall(r"<Iteration>(.*?)</Iteration>", xml_text, re.DOTALL)
        if not iteration_blocks:
            iteration_blocks = [xml_text]

        for iter_block in iteration_blocks:
            q_id_match = re.search(r"<Iteration_query-def>(.*?)</Iteration_query-def>", iter_block)
            if not q_id_match:
                q_id_match = re.search(r"<Iteration_query-ID>(.*?)</Iteration_query-ID>", iter_block)
            query_id = q_id_match.group(1).strip() if q_id_match else "query"

            q_len_match = re.search(r"<Iteration_query-len>(\d+)</Iteration_query-len>", iter_block)
            query_len = int(q_len_match.group(1)) if q_len_match else 0

            hit_blocks = re.findall(r"<Hit>(.*?)</Hit>", iter_block, re.DOTALL)
            for h_block in hit_blocks:
                s_id_match = re.search(r"<Hit_id>(.*?)</Hit_id>", h_block)
                s_def_match = re.search(r"<Hit_def>(.*?)</Hit_def>", h_block)
                s_len_match = re.search(r"<Hit_len>(\d+)</Hit_len>", h_block)

                subject_id = s_id_match.group(1).strip() if s_id_match else "subject"
                subject_title = s_def_match.group(1).strip() if s_def_match else ""
                subject_len = int(s_len_match.group(1)) if s_len_match else 0

                hsps = []
                hsp_blocks = re.findall(r"<Hsp>(.*?)</Hsp>", h_block, re.DOTALL)
                for hsp_b in hsp_blocks:
                    bit_score_m = re.search(r"<Hsp_bit-score>([\d.]+)</Hsp_bit-score>", hsp_b)
                    score_m = re.search(r"<Hsp_score>([\d.]+)</Hsp_score>", hsp_b)
                    evalue_m = re.search(r"<Hsp_evalue>([\d.eE+-]+)</Hsp_evalue>", hsp_b)
                    identity_m = re.search(r"<Hsp_identity>(\d+)</Hsp_identity>", hsp_b)
                    align_len_m = re.search(r"<Hsp_align-len>(\d+)</Hsp_align-len>", hsp_b)
                    gaps_m = re.search(r"<Hsp_gaps>(\d+)</Hsp_gaps>", hsp_b)
                    q_from_m = re.search(r"<Hsp_query-from>(\d+)</Hsp_query-from>", hsp_b)
                    q_to_m = re.search(r"<Hsp_query-to>(\d+)</Hsp_query-to>", hsp_b)
                    h_from_m = re.search(r"<Hsp_hit-from>(\d+)</Hsp_hit-from>", hsp_b)
                    h_to_m = re.search(r"<Hsp_hit-to>(\d+)</Hsp_hit-to>", hsp_b)

                    align_len = int(align_len_m.group(1)) if align_len_m else 0
                    ident_count = int(identity_m.group(1)) if identity_m else 0
                    ident_pct = (ident_count / align_len * 100.0) if align_len > 0 else 0.0

                    hsp = HSP(
                        score=float(score_m.group(1)) if score_m else 0.0,
                        bit_score=float(bit_score_m.group(1)) if bit_score_m else 0.0,
                        evalue=float(evalue_m.group(1)) if evalue_m else 1.0,
                        identity=round(ident_pct, 2),
                        alignment_length=align_len,
                        gaps=int(gaps_m.group(1)) if gaps_m else 0,
                        q_start=int(q_from_m.group(1)) if q_from_m else 1,
                        q_end=int(q_to_m.group(1)) if q_to_m else align_len,
                        s_start=int(h_from_m.group(1)) if h_from_m else 1,
                        s_end=int(h_to_m.group(1)) if h_to_m else align_len,
                    )
                    hsps.append(hsp)

                if hsps:
                    hits.append(
                        Hit(
                            query_id=query_id,
                            subject_id=subject_id,
                            subject_title=subject_title,
                            query_length=query_len,
                            subject_length=subject_len,
                            hsps=hsps,
                        )
                    )

        return hits


# ==============================================================================
# 3. KARLIN-ALTSCHUL STATISTICS & HIT FILTERING
# ==============================================================================

class KarlinAltschulStatistics:
    """Calculates Karlin-Altschul E-values and bit scores."""

    # Default Karlin-Altschul parameters for standard BLASTN/BLASTP scoring schemes
    # BLASTN (match=+1, mismatch=-2): lambda=1.37, K=0.711
    # BLASTP (BLOSUM62): lambda=0.267, K=0.041
    LAMBDA_BLASTP = 0.267
    K_BLASTP = 0.041

    @classmethod
    def raw_score_to_bit_score(cls, raw_score: float, lam: float = LAMBDA_BLASTP, k: float = K_BLASTP) -> float:
        """Convert raw alignment score S to bit score S' = (lambda * S - ln K) / ln 2."""
        return (lam * raw_score - math.log(k)) / math.log(2)

    @classmethod
    def calculate_evalue(
        cls,
        bit_score: float,
        query_len: int,
        db_total_letters: int,
    ) -> float:
        """Compute Expectation value: E = m * n * 2^(-S')."""
        search_space = max(query_len * db_total_letters, 1)
        exponent = -bit_score * math.log(2)
        return search_space * math.exp(exponent)


class HitFilter:
    """Filters BLAST hits according to statistical significance and alignment bounds."""

    @classmethod
    def filter_hits(
        cls,
        hits: List[Hit],
        criteria: Optional[FilterCriteria] = None,
    ) -> Dict[str, Any]:
        """Apply multi-attribute filtering and return retained hits with audit statistics."""
        crit = criteria or FilterCriteria()
        kept: List[Hit] = []
        stats = {
            "total_input_hits": len(hits),
            "retained_hits": 0,
            "rejected_evalue": 0,
            "rejected_bit_score": 0,
            "rejected_identity": 0,
            "rejected_query_coverage": 0,
            "rejected_alignment_length": 0,
            "rejected_gaps": 0,
        }

        for h in hits:
            best = h.best_hsp
            if not best:
                continue

            if best.evalue > crit.max_evalue:
                stats["rejected_evalue"] += 1
                continue
            if best.bit_score < crit.min_bit_score:
                stats["rejected_bit_score"] += 1
                continue
            if best.identity < crit.min_identity:
                stats["rejected_identity"] += 1
                continue
            if best.alignment_length < crit.min_alignment_length:
                stats["rejected_alignment_length"] += 1
                continue
            if crit.max_gaps is not None and best.gaps > crit.max_gaps:
                stats["rejected_gaps"] += 1
                continue
            if h.query_coverage < crit.min_query_coverage:
                stats["rejected_query_coverage"] += 1
                continue

            kept.append(h)

        stats["retained_hits"] = len(kept)
        return {
            "retained_hits": kept,
            "filter_statistics": stats,
            "retention_rate_pct": round(len(kept) / max(len(hits), 1) * 100.0, 2),
        }


# ==============================================================================
# 4. TAXONOMIC LOWEST COMMON ANCESTOR (LCA)
# ==============================================================================

class TaxonomicLCAEngine:
    """Bitscore-weighted Lowest Common Ancestor (LCA) classifier."""

    @classmethod
    def compute_lca(
        cls,
        hits: List[Hit],
        min_support_bitscore_fraction: float = 0.51,
    ) -> Dict[str, Any]:
        """
        Compute bitscore-weighted LCA across taxonomic lineages of hits for a query.
        """
        hits_with_lineage = [h for h in hits if h.taxonomic_lineage]
        if not hits_with_lineage:
            return {
                "lca_taxon": None,
                "lca_rank_level": 0,
                "lineage_to_lca": [],
                "support_fraction": 0.0,
                "hits_evaluated": len(hits),
            }

        total_bit_score = sum(h.best_hsp.bit_score if h.best_hsp else 1.0 for h in hits_with_lineage)
        max_depth = max(len(h.taxonomic_lineage) for h in hits_with_lineage)

        assigned_lineage: List[str] = []
        lca_support = 0.0

        for depth in range(max_depth):
            rank_votes: Dict[str, float] = {}
            for h in hits_with_lineage:
                if depth < len(h.taxonomic_lineage):
                    taxon = h.taxonomic_lineage[depth]
                    score = h.best_hsp.bit_score if h.best_hsp else 1.0
                    rank_votes[taxon] = rank_votes.get(taxon, 0.0) + score

            best_taxon = None
            best_score = 0.0
            for taxon, score in rank_votes.items():
                if score > best_score:
                    best_score = score
                    best_taxon = taxon

            frac = best_score / max(total_bit_score, 1e-9)
            if best_taxon and frac >= min_support_bitscore_fraction:
                assigned_lineage.append(best_taxon)
                lca_support = frac
            else:
                break

        lca_taxon = assigned_lineage[-1] if assigned_lineage else None

        return {
            "lca_taxon": lca_taxon,
            "lca_rank_level": len(assigned_lineage),
            "lineage_to_lca": assigned_lineage,
            "support_fraction": round(lca_support, 4),
            "hits_evaluated": len(hits_with_lineage),
        }

    @classmethod
    def export_newick_tree(cls, query_lca_map: Dict[str, Dict[str, Any]]) -> str:
        """Convert a map of {query_id: lca_result} into a Newick taxonomic tree string."""
        tree: Dict[str, Any] = {}

        for q_id, res in query_lca_map.items():
            lineage = res.get("lineage_to_lca", [])
            if not lineage:
                continue
            curr = tree
            for taxon in lineage:
                curr = curr.setdefault(taxon, {})
            # Attach query leaf
            curr[f"'{q_id}'"] = {}

        def _to_newick(node: Dict[str, Any]) -> str:
            if not node:
                return ""
            children = []
            for name, subtree in node.items():
                child_str = _to_newick(subtree)
                if child_str:
                    children.append(f"({child_str}){name}")
                else:
                    children.append(name)
            return ",".join(children)

        content = _to_newick(tree)
        return f"({content});" if content else "();"


# ==============================================================================
# 5. SIX-FRAME TRANSLATION & ORF DETECTOR
# ==============================================================================

class SequenceTranslator:
    """Six-frame DNA sequence translator and ORF extractor."""

    @staticmethod
    def reverse_complement(dna_seq: str) -> str:
        """Compute reverse complement of nucleotide sequence."""
        return "".join(DNA_COMPLEMENT.get(base, "N") for base in reversed(dna_seq.strip()))

    @classmethod
    def translate_frame(cls, dna_seq: str, frame: int) -> str:
        """
        Translate nucleotide sequence in specified frame (1..6).
        Frames 1-3: Forward 5'->3' offset by (frame - 1).
        Frames 4-6: Reverse complement offset by (frame - 4).
        """
        if frame not in range(1, 7):
            raise ValueError("Frame must be an integer between 1 and 6")

        seq = dna_seq.upper().replace(" ", "").replace("\n", "")
        if frame > 3:
            seq = cls.reverse_complement(seq)
            offset = frame - 4
        else:
            offset = frame - 1

        peptides = []
        for i in range(offset, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            aa = STANDARD_CODON_TABLE.get(codon, "X")
            peptides.append(aa)

        return "".join(peptides)

    @classmethod
    def translate_six_frames(cls, dna_seq: str) -> Dict[int, str]:
        """Translate DNA sequence across all 6 reading frames."""
        return {f: cls.translate_frame(dna_seq, f) for f in range(1, 7)}

    @classmethod
    def find_orfs(
        cls,
        dna_seq: str,
        min_length_aa: int = 30,
        require_start_codon: bool = True,
    ) -> List[Dict[str, Any]]:
        """
        Identify open reading frames (ORFs) across all 6 frames.
        """
        orfs: List[Dict[str, Any]] = []
        for frame in range(1, 7):
            protein = cls.translate_frame(dna_seq, frame)
            idx = 0
            while idx < len(protein):
                if require_start_codon:
                    if protein[idx] == "M":
                        stop_idx = protein.find("*", idx)
                        if stop_idx == -1:
                            stop_idx = len(protein)
                        orf_seq = protein[idx:stop_idx]
                        if len(orf_seq) >= min_length_aa:
                            orfs.append({
                                "frame": frame,
                                "start_aa": idx + 1,
                                "end_aa": stop_idx,
                                "length_aa": len(orf_seq),
                                "sequence": orf_seq,
                            })
                        idx = stop_idx + 1
                    else:
                        idx += 1
                else:
                    stop_idx = protein.find("*", idx)
                    if stop_idx == -1:
                        stop_idx = len(protein)
                    orf_seq = protein[idx:stop_idx]
                    if len(orf_seq) >= min_length_aa:
                        orfs.append({
                            "frame": frame,
                            "start_aa": idx + 1,
                            "end_aa": stop_idx,
                            "length_aa": len(orf_seq),
                            "sequence": orf_seq,
                        })
                    idx = stop_idx + 1

        orfs.sort(key=lambda x: x["length_aa"], reverse=True)
        return orfs
