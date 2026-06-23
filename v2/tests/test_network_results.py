"""Tests for the miR-Network graph shaping (app_v1.network_results).

Focused on the PURE graph logic -- compute_network() and strip_lncrna_version()
-- which need no prediction output. The ceRNA "co-target bridge" semantics: a
(gene, lncRNA) pair survives only if some miRNA is predicted to target BOTH.

Run with:
    python3 -m unittest v2.tests.test_network_results
"""

import os
import sys
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app_v1.network_results import (  # noqa: E402
    compute_network,
    strip_lncrna_version,
)


class StripLncrnaVersionTests(unittest.TestCase):
    def test_strips_ensembl_version(self):
        self.assertEqual(strip_lncrna_version("ENST00000456328.2"), "ENST00000456328")

    def test_strips_flybase_version(self):
        self.assertEqual(strip_lncrna_version("FBtr0347114.1"), "FBtr0347114")

    def test_keeps_wormbase_isoform_suffix(self):
        # WormBase '.27' is part of the name, not a version -- must survive.
        self.assertEqual(strip_lncrna_version("Y51H4A.27"), "Y51H4A.27")

    def test_none_passthrough(self):
        self.assertIsNone(strip_lncrna_version(None))


class ComputeNetworkPairsModeTests(unittest.TestCase):
    """Pairs mode: only user (gene, lncRNA) pairs with a shared miRNA survive."""

    def setUp(self):
        # gene NM_1 hit by miR-A (2 tools), miR-B (1 tool); NM_2 hit by miR-C.
        self.gene_edges = {
            "NM_1": {"miR-A": {"miRanda", "PITA"}, "miR-B": {"DMISO"}},
            "NM_2": {"miR-C": {"miRanda"}},
        }
        # lncRNA ENST1 hit by miR-A, miR-C; ENST2 hit by miR-Z only.
        self.lncrna_edges = {
            "ENST1": {"miR-A": {"RNAhybrid"}, "miR-C": {"PITA"}},
            "ENST2": {"miR-Z": {"DMISO"}},
        }
        self.labels = {("NM_1", "gene"): ("TP53", "tumor protein p53"),
                       ("NM_2", "gene"): ("MYC", "MYC proto-oncogene")}

    def test_bridge_pair_kept_with_shared_mirna(self):
        pairs = [{"gene": "TP53", "gene_refseqs": ["NM_1"], "lncrna": "ENST1"}]
        net = compute_network(self.gene_edges, self.lncrna_edges,
                              pairs=pairs, labels=self.labels)
        self.assertEqual(net["mode"], "pairs")
        self.assertEqual(len(net["pairs"]), 1)
        kept = net["pairs"][0]
        # miR-A bridges NM_1 and ENST1; miR-B/miR-C do not target both.
        self.assertEqual(kept["bridge_mirnas"], ["miR-A"])
        self.assertEqual(kept["gene_label"], "TP53")
        node_ids = {n["id"] for n in net["nodes"]}
        self.assertEqual(node_ids, {"NM_1", "miR-A", "ENST1"})

    def test_pair_without_shared_mirna_dropped(self):
        # NM_1's miRNAs {A,B} vs ENST2's {Z}: no overlap -> pair dropped, empty graph.
        pairs = [{"gene": "TP53", "gene_refseqs": ["NM_1"], "lncrna": "ENST2"}]
        net = compute_network(self.gene_edges, self.lncrna_edges,
                              pairs=pairs, labels=self.labels)
        self.assertEqual(net["pairs"], [])
        self.assertEqual(net["nodes"], [])
        self.assertEqual(net["edges"], [])

    def test_edges_carry_full_tool_lists(self):
        pairs = [{"gene": "TP53", "gene_refseqs": ["NM_1"], "lncrna": "ENST1"}]
        net = compute_network(self.gene_edges, self.lncrna_edges,
                              pairs=pairs, labels=self.labels)
        gene_edge = [e for e in net["edges"] if e["side"] == "gene"][0]
        self.assertEqual(gene_edge["source"], "NM_1")
        self.assertEqual(gene_edge["target"], "miR-A")
        # tools NM_1<-miR-A = {miRanda, PITA}; client thresholds on tool_count.
        self.assertEqual(gene_edge["tools"], ["PITA", "miRanda"])
        self.assertEqual(gene_edge["tool_count"], 2)
        lncrna_edge = [e for e in net["edges"] if e["side"] == "lncrna"][0]
        self.assertEqual(lncrna_edge["tool_count"], 1)

    def test_symbol_expands_to_multiple_refseqs(self):
        # A gene symbol resolving to NM_1 and NM_2: both checked against the lncRNA.
        gene_edges = {
            "NM_1": {"miR-A": {"miRanda"}},
            "NM_2": {"miR-A": {"PITA"}},
        }
        lncrna_edges = {"ENST1": {"miR-A": {"DMISO"}}}
        pairs = [{"gene": "FAM", "gene_refseqs": ["NM_1", "NM_2"], "lncrna": "ENST1"}]
        net = compute_network(gene_edges, lncrna_edges, pairs=pairs)
        gene_node_ids = {n["id"] for n in net["nodes"] if n["type"] == "gene"}
        self.assertEqual(gene_node_ids, {"NM_1", "NM_2"})
        self.assertEqual(len(net["pairs"]), 2)

    def test_lncrna_version_normalized_against_parsed_ids(self):
        # User supplies a versioned lncRNA; parsed edges are unversioned.
        pairs = [{"gene": "TP53", "gene_refseqs": ["NM_1"], "lncrna": "ENST1.4"}]
        net = compute_network(self.gene_edges, self.lncrna_edges,
                              pairs=pairs, labels=self.labels)
        self.assertEqual(len(net["pairs"]), 1)
        self.assertIn("ENST1", {n["id"] for n in net["nodes"]})


class ComputeNetworkDiscoveryModeTests(unittest.TestCase):
    """No pairs -> bounded discovery over top-degree genes/lncRNAs."""

    def test_discovery_builds_bridges_among_top_targets(self):
        gene_edges = {"NM_1": {"miR-A": {"miRanda"}}}
        lncrna_edges = {"ENST1": {"miR-A": {"PITA"}}}
        net = compute_network(gene_edges, lncrna_edges, pairs=None)
        self.assertEqual(net["mode"], "discovery")
        self.assertEqual(net["pairs"][0]["bridge_mirnas"], ["miR-A"])
        self.assertEqual({n["id"] for n in net["nodes"]}, {"NM_1", "miR-A", "ENST1"})

    def test_discovery_respects_top_n_caps(self):
        # Two genes; only the highest-degree one is kept when top_genes=1.
        gene_edges = {
            "NM_hi": {"miR-A": {"x"}, "miR-B": {"x"}},  # degree 2
            "NM_lo": {"miR-A": {"x"}},                    # degree 1
        }
        lncrna_edges = {"ENST1": {"miR-A": {"y"}, "miR-B": {"y"}}}
        net = compute_network(gene_edges, lncrna_edges, pairs=None, top_genes=1)
        gene_node_ids = {n["id"] for n in net["nodes"] if n["type"] == "gene"}
        self.assertEqual(gene_node_ids, {"NM_hi"})

    def test_empty_inputs_yield_empty_graph(self):
        net = compute_network({}, {}, pairs=None)
        self.assertEqual(net["nodes"], [])
        self.assertEqual(net["edges"], [])
        self.assertEqual(net["summary"]["edge_count"], 0)


if __name__ == "__main__":
    unittest.main()
