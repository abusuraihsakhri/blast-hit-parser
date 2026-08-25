"""Blast XML Parser: parse BLAST XML output, extract HSPs, compute E-value statistics."""
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
import math
import re


@dataclass
class HighScoringPair:
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    score: float
    e_value: float
    identity: float
    alignment_length: int
    gaps: int = 0
    strand: str = "+"


@dataclass
class BlastHit:
    query_id: str
    subject_id: str
    subject_description: str
    hsps: List[HighScoringPair]
    total_score: float = 0.0
    query_coverage: float = 0.0
    percent_identity: float = 0.0


class BlastXMLParser:
    """Parse BLAST XML output (text-based, no lxml dependency)."""

    def parse_xml_text(self, xml_text: str) -> Dict[str, Any]:
        """Parse BLAST XML output from text string."""
        hits = []
        hit_blocks = re.findall(r'<Hit>.*?</Hit>', xml_text, re.DOTALL)

        for block in hit_blocks:
            hit = self._parse_hit_block(block)
            if hit:
                hits.append(hit)

        return {
            "hits": hits,
            "num_hits": len(hits),
            "query_ids": list(set(h.query_id for h in hits)),
        }

    def _parse_hit_block(self, block: str) -> Optional[BlastHit]:
        query_id_m = re.search(r'<Hit_id>(.*?)</Hit_id>', block)
        desc_m = re.search(r'<Hit_def>(.*?)</Hit_def>', block)
        if not query_id_m:
            return None

        query_id = query_id_m.group(1).strip()
        subject_desc = desc_m.group(1).strip() if desc_m else ""

        hsp_blocks = re.findall(r'<Hsp>.*?</Hsp>', block, re.DOTALL)
        hsps = []
        for hb in hsp_blocks:
            hsp = self._parse_hsp(hb)
            if hsp:
                hsps.append(hsp)

        if not hsps:
            return None

        total_score = sum(h.score for h in hsps)
        subject_id = query_id

        max_end = max((h.query_end for h in hsps), default=0)
        min_start = min((h.query_start for h in hsps), default=0)
        query_len = max_end - min_start + 1 if max_end > min_start else 1
        query_coverage = sum(h.alignment_length for h in hsps) / max(query_len, 1)

        total_identity = sum(h.identity * h.alignment_length for h in hsps)
        total_len = sum(h.alignment_length for h in hsps)
        percent_identity = total_identity / max(total_len, 1)

        return BlastHit(
            query_id=query_id, subject_id=subject_id,
            subject_description=subject_desc, hsps=hsps,
            total_score=total_score, query_coverage=min(query_coverage, 1.0),
            percent_identity=round(percent_identity * 100, 2),
        )

    def _parse_hsp(self, block: str) -> Optional[HighScoringPair]:
        def extract_float(tag):
            m = re.search(f'<{tag}>(.*?)</{tag}>', block)
            return float(m.group(1)) if m else 0.0

        def extract_int(tag):
            m = re.search(f'<{tag}>(.*?)</{tag}>', block)
            return int(m.group(1)) if m else 0

        score = extract_float('Hsp_score')
        e_value = extract_float('Hsp_evalue')
        identity = extract_int('Hsp_identity')
        align_len = extract_int('Hsp_align-len')
        gaps = extract_int('Hsp_gaps')

        q_from = extract_int('Hsp_query-from')
        q_to = extract_int('Hsp_query-to')
        h_from = extract_int('Hsp_hit-from')
        h_to = extract_int('Hsp_hit-to')

        if align_len == 0:
            return None

        return HighScoringPair(
            query_start=q_from, query_end=q_to,
            subject_start=h_from, subject_end=h_to,
            score=score, e_value=e_value,
            identity=identity / align_len if align_len > 0 else 0,
            alignment_length=align_len, gaps=gaps,
        )

    def evalue_statistics(self, hits: List[BlastHit]) -> Dict[str, Any]:
        """Compute E-value distribution statistics."""
        all_evalues = []
        for h in hits:
            for hsp in h.hsps:
                if hsp.e_value > 0:
                    all_evalues.append(hsp.e_value)

        if not all_evalues:
            return {"status": "no_hits"}

        log_evals = [math.log10(e) for e in all_evalues]
        mean_log = sum(log_evals) / len(log_evals)
        std_log = math.sqrt(sum((x - mean_log) ** 2 for x in log_evals) / max(len(log_evals) - 1, 1))

        significant = [e for e in all_evalues if e < 1e-5]

        return {
            "total_hsps": len(all_evalues),
            "min_evalue": min(all_evalues),
            "max_evalue": max(all_evalues),
            "mean_log10_evalue": round(mean_log, 4),
            "std_log10_evalue": round(std_log, 4),
            "significant_hits_count": len(significant),
            "significant_percent": round(len(significant) / len(all_evalues) * 100, 1),
        }
