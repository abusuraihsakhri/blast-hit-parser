"""Blast Translation: translate nucleotide queries, reverse complement, codon table support."""
from typing import Dict, Any, Optional


STANDARD_CODON_TABLE = {
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

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


class SequenceTranslator:
    """Translate nucleotide sequences in all six reading frames."""

    def reverse_complement(self, seq: str) -> str:
        return "".join(COMPLEMENT.get(c.upper(), "N") for c in reversed(seq))

    def translate_frame(self, seq: str, frame: int) -> str:
        """Translate in one reading frame (1-6: 1-3 forward, 4-6 reverse)."""
        if frame < 1 or frame > 6:
            raise ValueError("Frame must be 1-6")

        if frame > 3:
            seq = self.reverse_complement(seq)
            frame -= 3

        start = frame - 1
        protein = []
        for i in range(start, len(seq) - 2, 3):
            codon = seq[i:i+3].upper()
            aa = STANDARD_CODON_TABLE.get(codon, "X")
            protein.append(aa)
        return "".join(protein)

    def translate_six_frames(self, seq: str) -> Dict[str, str]:
        """Translate all six reading frames."""
        return {f"frame_{f}": self.translate_frame(seq, f) for f in range(1, 7)}

    def find_open_reading_frames(self, seq: str, min_length_aa: int = 30) -> Dict[str, Any]:
        """Find all ORFs above minimum length."""
        orfs = []
        for frame in range(1, 7):
            protein = self.translate_frame(seq, frame)
            orfs.extend(self._extract_orfs(protein, frame, min_length_aa))

        orfs.sort(key=lambda x: -x["length_aa"])
        return {"orfs": orfs, "num_orfs": len(orfs)}

    def _extract_orfs(self, protein: str, frame: int, min_length: int) -> list:
        orfs = []
        i = 0
        while i < len(protein):
            if protein[i] == "M":
                j = i
                while j < len(protein) and protein[j] != "*":
                    j += 1
                length = j - i
                if length >= min_length:
                    orfs.append({
                        "frame": frame,
                        "start_aa": i + 1,
                        "end_aa": j,
                        "length_aa": length,
                        "sequence": protein[i:j],
                        "ends_with_stop": j < len(protein) and protein[j] == "*",
                    })
                i = j + 1
            else:
                i += 1
        return orfs

    def gc_content(self, seq: str) -> Dict[str, Any]:
        """Compute GC content and sliding window GC."""
        seq_upper = seq.upper()
        gc = sum(1 for c in seq_upper if c in "GC")
        total = len(seq_upper) if seq_upper else 1
        overall_gc = gc / total

        window_size = 100
        window_gc = []
        for i in range(0, len(seq_upper) - window_size + 1, window_size // 2):
            w = seq_upper[i:i + window_size]
            wgc = sum(1 for c in w if c in "GC") / max(len(w), 1)
            window_gc.append(round(wgc, 4))

        return {
            "overall_gc": round(overall_gc, 4),
            "length": len(seq_upper),
            "gc_windows": window_gc,
            "gc_skew": round((sum(1 for c in seq_upper if c == "G") - sum(1 for c in seq_upper if c == "C")) / max(total, 1), 4),
        }
