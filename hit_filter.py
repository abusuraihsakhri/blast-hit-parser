"""Blast Hit Filtering and Clustering: threshold filters, query-subject clustering, overlap resolution."""
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field


@dataclass
class FilterCriteria:
    max_evalue: float = 1e-5
    min_identity_percent: float = 30.0
    min_alignment_length: int = 50
    min_score: float = 40.0
    max_gaps: int = 10


class HitFilter:
    """Filter BLAST hits by configurable criteria."""

    def filter_hits(self, hits: List[Dict], criteria: FilterCriteria) -> Dict[str, Any]:
        """Apply multi-criteria filtering."""
        filtered = []
        rejected = {"evalue": 0, "identity": 0, "alignment_length": 0, "score": 0, "gaps": 0}

        for hit in hits:
            hsps = hit.get("hsps", [])
            passing_hsps = []
            for hsp in hsps:
                reasons = []
                if hsp.get("e_value", 1) > criteria.max_evalue:
                    reasons.append("evalue")
                if hsp.get("identity", 0) * 100 < criteria.min_identity_percent:
                    reasons.append("identity")
                if hsp.get("alignment_length", 0) < criteria.min_alignment_length:
                    reasons.append("alignment_length")
                if hsp.get("score", 0) < criteria.min_score:
                    reasons.append("score")
                if hsp.get("gaps", 0) > criteria.max_gaps:
                    reasons.append("gaps")

                if not reasons:
                    passing_hsps.append(hsp)
                else:
                    for r in reasons:
                        rejected[r] += 1

            if passing_hsps:
                filtered.append({**hit, "hsps": passing_hsps})

        return {
            "total_hits": len(hits),
            "passing_hits": len(filtered),
            "rejected_counts": rejected,
            "filtered_hits": filtered,
        }


class HitClusterer:
    """Cluster BLAST hits by query/subject overlap."""

    def cluster_by_subject(self, hits: List[Dict]) -> Dict[str, Any]:
        """Group hits by subject accession."""
        clusters: Dict[str, List[Dict]] = {}
        for hit in hits:
            subj = hit.get("subject_id", "unknown")
            clusters.setdefault(subj, []).append(hit)

        result = []
        for subj, cluster_hits in sorted(clusters.items(), key=lambda x: -len(x[1])):
            best_score = max(h.get("total_score", 0) for h in cluster_hits)
            result.append({
                "subject_id": subj,
                "num_hits": len(cluster_hits),
                "best_score": best_score,
                "query_ids": list(set(h.get("query_id", "") for h in cluster_hits)),
            })

        return {"clusters": result, "num_clusters": len(result)}

    def resolve_overlaps(self, hits: List[Dict]) -> List[Dict]:
        """Remove overlapping HSPs, keeping highest-scoring per query region."""
        sorted_hits = sorted(hits, key=lambda h: -h.get("total_score", 0))
        covered_regions: List[tuple] = []
        resolved = []

        for hit in sorted_hits:
            hsps = hit.get("hsps", [])
            non_overlapping = []
            for hsp in hsps:
                qs, qe = hsp.get("query_start", 0), hsp.get("query_end", 0)
                region = (min(qs, qe), max(qs, qe))
                overlaps = False
                for cs, ce in covered_regions:
                    if region[0] < ce and region[1] > cs:
                        overlaps = True
                        break
                if not overlaps:
                    non_overlapping.append(hsp)
                    covered_regions.append(region)

            if non_overlapping:
                resolved.append({**hit, "hsps": non_overlapping})

        return resolved
