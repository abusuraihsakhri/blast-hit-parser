"""
Unit Test Suite for BLAST Hit Parser & Taxonomic LCA Engine
============================================================
Comprehensive test suite covering:
  - BLAST tabular (outfmt 6 / m8) parsing with custom & standard headers
  - BLAST XML (outfmt 5) parsing
  - Karlin-Altschul E-value and bit score calculations
  - Multi-criteria hit filtering and query coverage interval merging
  - Bitscore-weighted taxonomic Lowest Common Ancestor (LCA) assignment
  - Newick tree serialization
  - Six-frame translation and ORF detection
  - CLI argument handling and JSON output
"""

import io
import json
import unittest
from contextlib import redirect_stdout

from blast_hit_parser import (
    BlastTabularParser,
    BlastXMLParser,
    FilterCriteria,
    HSP,
    Hit,
    HitFilter,
    KarlinAltschulStatistics,
    SequenceTranslator,
    TaxonomicLCAEngine,
)
from cli import main


class TestBlastTabularParser(unittest.TestCase):
    """Test BLAST tabular format parsing."""

    SAMPLE_TABULAR = """
# BLASTN 2.14.0+
# Query: seq1_16S_rRNA
# Fields: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
seq1\tNR_024570.1_Ecoli\t99.50\t1500\t7\t1\t1\t1500\t1\t1500\t0.0\t2750.0
seq1\tNR_042817.1_Shigella\t99.20\t1500\t12\t0\t1\t1500\t1\t1500\t0.0\t2720.0
seq1\tNR_074910.1_Salmonella\t94.80\t1480\t75\t2\t10\t1485\t5\t1480\t1e-150\t1850.0
seq2\tNR_112001.1_Staph\t98.00\t500\t10\t0\t1\t500\t1\t500\t0.0\t920.0
"""

    def test_parse_tabular_basic(self):
        hits = BlastTabularParser.parse_text(self.SAMPLE_TABULAR)
        self.assertEqual(len(hits), 4)
        h1 = hits[0]
        self.assertEqual(h1.query_id, "seq1")
        self.assertEqual(h1.subject_id, "NR_024570.1_Ecoli")
        self.assertEqual(h1.best_hsp.identity, 99.50)
        self.assertEqual(h1.best_hsp.alignment_length, 1500)
        self.assertEqual(h1.best_hsp.evalue, 0.0)

    def test_multi_hsp_grouping(self):
        multi_hsp_text = """
seq1\tsubjA\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-40\t180.0
seq1\tsubjA\t90.0\t100\t10\t0\t150\t250\t101\t200\t1e-30\t140.0
"""
        hits = BlastTabularParser.parse_text(multi_hsp_text)
        self.assertEqual(len(hits), 1)
        self.assertEqual(len(hits[0].hsps), 2)
        self.assertEqual(hits[0].total_bit_score, 320.0)

    def test_query_coverage_interval_merging(self):
        # Query length 300; two HSPs covering 1-100 and 80-200 -> total covered 1-200 = 200 bp -> cov = 200/300 = 66.67%
        hsp1 = HSP(100, 100, 1e-20, 95.0, 100, q_start=1, q_end=100)
        hsp2 = HSP(100, 100, 1e-20, 95.0, 120, q_start=80, q_end=200)
        hit = Hit(query_id="q1", subject_id="s1", query_length=300, hsps=[hsp1, hsp2])
        self.assertAlmostEqual(hit.query_coverage, 200.0 / 300.0, places=4)

    def test_custom_taxonomy_lineage_parsing(self):
        tax_text = "q1\ts1\t99.0\t100\t1\t0\t1\t100\t1\t100\t0.0\t200.0\tBacteria;Proteobacteria;Escherichia"
        fields = BlastTabularParser.STANDARD_FIELDS + ["lineage"]
        hits = BlastTabularParser.parse_text(tax_text, custom_fields=fields)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].taxonomic_lineage, ["Bacteria", "Proteobacteria", "Escherichia"])


class TestBlastXMLParser(unittest.TestCase):
    """Test BLAST XML output parser."""

    SAMPLE_XML = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>test_gene_X</Iteration_query-def>
      <Iteration_query-len>600</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>ref|NP_001234.1|</Hit_id>
          <Hit_def>DNA polymerase [Escherichia coli]</Hit_def>
          <Hit_accession>NP_001234</Hit_accession>
          <Hit_len>580</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>450.5</Hsp_bit-score>
              <Hsp_score>1158</Hsp_score>
              <Hsp_evalue>2.5e-120</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>580</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>580</Hsp_hit-to>
              <Hsp_identity>560</Hsp_identity>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>580</Hsp_align-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""

    def test_xml_parsing(self):
        hits = BlastXMLParser.parse_xml_text(self.SAMPLE_XML)
        self.assertEqual(len(hits), 1)
        h = hits[0]
        self.assertEqual(h.query_id, "test_gene_X")
        self.assertEqual(h.subject_id, "ref|NP_001234.1|")
        self.assertEqual(h.query_length, 600)
        self.assertEqual(len(h.hsps), 1)
        hsp = h.hsps[0]
        self.assertEqual(hsp.bit_score, 450.5)
        self.assertEqual(hsp.evalue, 2.5e-120)
        self.assertAlmostEqual(hsp.identity, (560 / 580) * 100.0, places=2)
        self.assertEqual(hsp.gaps, 0)


class TestKarlinAltschulStatistics(unittest.TestCase):
    """Test Karlin-Altschul statistical conversions."""

    def test_score_to_bit_score(self):
        raw_score = 100.0
        bit_score = KarlinAltschulStatistics.raw_score_to_bit_score(raw_score)
        self.assertGreater(bit_score, 0.0)

    def test_evalue_calculation_scaling(self):
        ev1 = KarlinAltschulStatistics.calculate_evalue(bit_score=50.0, query_len=300, db_total_letters=1000000)
        ev2 = KarlinAltschulStatistics.calculate_evalue(bit_score=100.0, query_len=300, db_total_letters=1000000)
        # Higher bit score should produce drastically lower E-value
        self.assertGreater(ev1, ev2)


class TestHitFilter(unittest.TestCase):
    """Test filtering of BLAST hits."""

    def setUp(self):
        h1 = Hit(
            "q1", "s_pass", query_length=100,
            hsps=[HSP(score=200, bit_score=200, evalue=1e-50, identity=95.0, alignment_length=95, q_start=1, q_end=95)]
        )
        h2 = Hit(
            "q1", "s_low_ident", query_length=100,
            hsps=[HSP(score=100, bit_score=80, evalue=1e-20, identity=50.0, alignment_length=95, q_start=1, q_end=95)]
        )
        h3 = Hit(
            "q1", "s_high_evalue", query_length=100,
            hsps=[HSP(score=40, bit_score=30, evalue=0.05, identity=90.0, alignment_length=95, q_start=1, q_end=95)]
        )
        self.hits = [h1, h2, h3]

    def test_filter_hits_criteria(self):
        crit = FilterCriteria(max_evalue=1e-5, min_bit_score=50.0, min_identity=70.0)
        res = HitFilter.filter_hits(self.hits, crit)
        self.assertEqual(len(res["retained_hits"]), 1)
        self.assertEqual(res["retained_hits"][0].subject_id, "s_pass")
        self.assertEqual(res["filter_statistics"]["rejected_identity"], 1)
        self.assertEqual(res["filter_statistics"]["rejected_evalue"], 1)


class TestTaxonomicLCA(unittest.TestCase):
    """Test Lowest Common Ancestor assignment and Newick tree export."""

    def test_unanimous_lca(self):
        h1 = Hit("q1", "s1", hsps=[HSP(100, 100, 0.0, 99.0, 100)],
                 taxonomic_lineage=["Bacteria", "Proteobacteria", "Enterobacterales", "Escherichia", "E_coli"])
        h2 = Hit("q1", "s2", hsps=[HSP(90, 90, 0.0, 98.0, 100)],
                 taxonomic_lineage=["Bacteria", "Proteobacteria", "Enterobacterales", "Escherichia", "E_coli"])
        lca = TaxonomicLCAEngine.compute_lca([h1, h2])
        self.assertEqual(lca["lca_taxon"], "E_coli")
        self.assertEqual(lca["lca_rank_level"], 5)
        self.assertEqual(lca["support_fraction"], 1.0)

    def test_divergent_genus_lca(self):
        h1 = Hit("q1", "s1", hsps=[HSP(100, 100, 0.0, 99.0, 100)],
                 taxonomic_lineage=["Bacteria", "Proteobacteria", "Enterobacterales", "Escherichia"])
        h2 = Hit("q1", "s2", hsps=[HSP(100, 100, 0.0, 99.0, 100)],
                 taxonomic_lineage=["Bacteria", "Proteobacteria", "Enterobacterales", "Salmonella"])
        lca = TaxonomicLCAEngine.compute_lca([h1, h2])
        # Diverges at genus level -> LCA is Enterobacterales
        self.assertEqual(lca["lca_taxon"], "Enterobacterales")
        self.assertEqual(lca["lineage_to_lca"], ["Bacteria", "Proteobacteria", "Enterobacterales"])

    def test_empty_lineage_lca(self):
        h1 = Hit("q1", "s1", hsps=[HSP(100, 100, 0.0, 99.0, 100)], taxonomic_lineage=[])
        lca = TaxonomicLCAEngine.compute_lca([h1])
        self.assertIsNone(lca["lca_taxon"])

    def test_newick_tree_generation(self):
        query_map = {
            "q1": {"lineage_to_lca": ["Bacteria", "Firmicutes", "Bacillus"]},
            "q2": {"lineage_to_lca": ["Bacteria", "Firmicutes", "Staphylococcus"]},
        }
        newick = TaxonomicLCAEngine.export_newick_tree(query_map)
        self.assertTrue(newick.endswith(";"))
        self.assertIn("Bacteria", newick)
        self.assertIn("Bacillus", newick)


class TestSequenceTranslator(unittest.TestCase):
    """Test 6-frame translation and ORF finder."""

    def test_forward_frame_1(self):
        # ATG (M) GCT (A) GGT (G) TAA (*)
        dna = "ATGGCTGGGTAA"
        prot = SequenceTranslator.translate_frame(dna, 1)
        self.assertEqual(prot, "MAG*")

    def test_reverse_complement_frame(self):
        # Forward: TTACCCAGCCAT -> Reverse complement: ATGGCTGGGTAA -> Frame 4: MAG*
        dna = "TTACCCAGCCAT"
        prot = SequenceTranslator.translate_frame(dna, 4)
        self.assertEqual(prot, "MAG*")

    def test_find_orfs(self):
        # Sequence containing an ORF of length 5: ATG GCT GCT GCT TAA
        dna = "CCCATGGCTGCTGCTTAAGGG"
        orfs = SequenceTranslator.find_orfs(dna, min_length_aa=4, require_start_codon=True)
        self.assertEqual(len(orfs), 1)
        self.assertEqual(orfs[0]["sequence"], "MAAA")

    def test_invalid_frame_error(self):
        with self.assertRaises(ValueError):
            SequenceTranslator.translate_frame("ATG", 7)


class TestCLIExecution(unittest.TestCase):
    """Test CLI commands and JSON format."""

    SAMPLE_TABULAR = "q1\ts1\t99.0\t100\t1\t0\t1\t100\t1\t100\t1e-20\t200.0\n"

    def test_cli_parse_tabular_json(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            ret = main(["parse-tabular", "--text", self.SAMPLE_TABULAR, "--json"])
        self.assertEqual(ret, 0)
        data = json.loads(buf.getvalue())
        self.assertEqual(len(data), 1)
        self.assertEqual(data[0]["query_id"], "q1")

    def test_cli_filter_json(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            ret = main(["filter", "--text", self.SAMPLE_TABULAR, "--max-evalue", "1e-5", "--json"])
        self.assertEqual(ret, 0)
        data = json.loads(buf.getvalue())
        self.assertEqual(data["retention_rate_pct"], 100.0)

    def test_cli_translate_json(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            ret = main(["translate", "--sequence", "ATGGCTTAA", "--json"])
        self.assertEqual(ret, 0)
        data = json.loads(buf.getvalue())
        self.assertIn("six_frames", data)

    def test_cli_lca_json(self):
        tax_text = "q1\ts1\t99.0\t100\t1\t0\t1\t100\t1\t100\t0.0\t200.0\tBacteria;Proteobacteria\n"
        buf = io.StringIO()
        with redirect_stdout(buf):
            ret = main(["lca", "--text", tax_text, "--json"])
        self.assertEqual(ret, 0)
        data = json.loads(buf.getvalue())
        self.assertIn("assignments", data)
        self.assertIn("newick_tree", data)

    def test_cli_parse_xml_json(self):
        xml_sample = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_query-def>q_test</Iteration_query-def>
      <Iteration_hits>
        <Hit>
          <Hit_id>s_test</Hit_id>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>150.0</Hsp_bit-score>
              <Hsp_score>300</Hsp_score>
              <Hsp_evalue>1e-35</Hsp_evalue>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_identity>98</Hsp_identity>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""
        buf = io.StringIO()
        with redirect_stdout(buf):
            ret = main(["parse-xml", "--text", xml_sample, "--json"])
        self.assertEqual(ret, 0)
        data = json.loads(buf.getvalue())
        self.assertEqual(len(data), 1)
        self.assertEqual(data[0]["query_id"], "q_test")


class TestEdgeCasesAndStatistics(unittest.TestCase):
    """Test mathematical and parsing edge cases."""

    def test_empty_tabular_input(self):
        hits = BlastTabularParser.parse_text("")
        self.assertEqual(len(hits), 0)

    def test_whitespace_delimited_tabular(self):
        line = "q1 s1 99.0 100 0 0 1 100 1 100 0.0 200.0"
        hits = BlastTabularParser.parse_text(line)
        self.assertEqual(len(hits), 1)

    def test_karlin_altschul_custom_lambda(self):
        raw = 50.0
        # BLASTN parameters: lambda=1.37, K=0.711
        bit = KarlinAltschulStatistics.raw_score_to_bit_score(raw, lam=1.37, k=0.711)
        self.assertGreater(bit, raw)

    def test_orf_without_start_codon_requirement(self):
        # A peptide sequence without Met at the start: GCT GCT GCT TAA -> AAA*
        dna = "GCTGCTGCTTAA"
        orfs = SequenceTranslator.find_orfs(dna, min_length_aa=3, require_start_codon=False)
        seqs = [o["sequence"] for o in orfs]
        self.assertIn("AAA", seqs)

    def test_filter_no_passing_hits(self):
        h = Hit("q1", "s1", hsps=[HSP(10, 10, 0.9, 30.0, 20)])
        crit = FilterCriteria(max_evalue=1e-10, min_bit_score=100.0)
        res = HitFilter.filter_hits([h], crit)
        self.assertEqual(len(res["retained_hits"]), 0)
        self.assertEqual(res["retention_rate_pct"], 0.0)

    def test_lca_low_support_truncation(self):
        h1 = Hit("q1", "s1", hsps=[HSP(100, 100, 0.0, 99.0, 100)],
                 taxonomic_lineage=["Bacteria", "Proteobacteria", "Escherichia"])
        h2 = Hit("q1", "s2", hsps=[HSP(200, 200, 0.0, 99.0, 100)],
                 taxonomic_lineage=["Bacteria", "Actinobacteria", "Mycobacterium"])
        lca = TaxonomicLCAEngine.compute_lca([h1, h2], min_support_bitscore_fraction=0.70)
        self.assertEqual(lca["lca_taxon"], "Bacteria")
        self.assertEqual(lca["lineage_to_lca"], ["Bacteria"])

    def test_reverse_complement_ambiguous_bases(self):
        # DNA containing N, U, and standard bases: A, T, G, C, U, N
        # In reverse: N, U, C, G, T, A -> Complements: N, A, G, C, A, T -> "NAGCAT"
        dna = "ATGCUN"
        rc = SequenceTranslator.reverse_complement(dna)
        self.assertEqual(rc, "NAGCAT")

    def test_query_coverage_zero_query_len(self):
        h = Hit("q1", "s1", query_length=0, hsps=[HSP(10, 10, 0.0, 90.0, 50)])
        self.assertEqual(h.query_coverage, 0.0)

    def test_hit_best_hsp_empty(self):
        h = Hit("q1", "s1", hsps=[])
        self.assertIsNone(h.best_hsp)
        self.assertEqual(h.best_evalue, float("inf"))


if __name__ == "__main__":
    unittest.main()
