"""Tests for per-sequence results in result.db / lncrna_result.db.

A job that adds a shift or a base modification runs every tool on the WT and
on each variant, and each run writes its own raw file named after the record's
header. result_db used to pour all of them into one gene_tools(gene_id, tool)
table, so the results table, the Venn/UpSet plots, the enrichment gene list and
the download showed the union: variant-only targets were unmarked, targets the
variant lost were hidden behind the persisting WT edge, and tool-agreement
counts were computed on sets no single sequence produced.

These tests stub the raw-file parsers and cover the aggregation layer only:
schema version 2, the WT default view, per-record and union views, the WT vs
variant comparison, the lazy rebuild of a version-1 database, and the CSV.

Run with:
    python3 -m unittest v2.tests.test_result_db_sequences
"""

import csv
import json
import os
import shutil
import sqlite3
import sys
import tempfile
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from app_v1 import result_db  # noqa: E402
from app_v1 import lncrna_results  # noqa: E402
from app_v1 import network_results  # noqa: E402
from app_v1.sequence_ids import (  # noqa: E402
    is_variant, mirna_node_id, mirna_node_label, sequence_record, variant_description,
)

WT = "hsa-miR-181d-5p"
VAR = "hsa-miR-181d-5p,6:U|A,modified"
WT_HEADER = "hsa-miR-181d-5p,WT"

# Sebastiano Di Bella's f5a73642 in miniature: with 6:U|A the variant loses
# the miRanda and PITA sites on MBNL1 (NM_021038) but keeps RNAhybrid, drops
# MYC altogether, keeps TP53 unchanged and gains EZH2.
PREDICTIONS = {
    WT_HEADER: {
        "miRanda":   ["NM_021038", "NM_000546"],
        "PITA":      ["NM_021038"],
        "RNAhybrid": ["NM_021038", "NM_000546", "NM_002467"],
    },
    VAR: {
        "miRanda":   ["NM_000546", "NM_004456"],
        "PITA":      [],
        "RNAhybrid": ["NM_021038", "NM_000546"],
    },
}
LABELS = {
    "NM_021038": ("MBNL1", "muscleblind like splicing regulator 1"),
    "NM_000546": ("TP53", "tumor protein p53"),
    "NM_002467": ("MYC", "MYC proto-oncogene"),
    "NM_004456": ("EZH2", "enhancer of zeste 2"),
}


def _write_params(output_dir, headers):
    seqs = []
    for h in headers:
        seqs.append({"header": h, "sequence": "AACAUUCAUUGUUGUCGGUGGGU",
                     "length": 23, "type": h.split(",")[-1]})
    with open(os.path.join(output_dir, "mirna_prediction_parameters.json"), "w") as fh:
        json.dump({"sequences": seqs, "tools": ["miRanda", "PITA", "RNAhybrid"],
                   "num_cores": 1}, fh)


def _fake_process_sequence(sequence, output_dir, **_kw):
    return {"header": sequence["header"],
            "prediction": PREDICTIONS.get(sequence["header"], {})}


def _tools(db_path, sequence=None):
    """{gene_id: set(tools)} straight from the table."""
    conn = sqlite3.connect(db_path)
    try:
        if sequence is None:
            rows = conn.execute("SELECT gene_id, tool FROM gene_tools").fetchall()
        else:
            rows = conn.execute("SELECT gene_id, tool FROM gene_tools WHERE sequence = ?",
                                (sequence,)).fetchall()
    finally:
        conn.close()
    out = {}
    for gene_id, tool in rows:
        out.setdefault(gene_id, set()).add(tool)
    return out


class _StubbedParsers(unittest.TestCase):
    """Stub everything in result_db that touches raw files or reference dbs."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="isotar_seq_")
        self.output_dir = os.path.join(self.tmp, "output")
        os.makedirs(self.output_dir)
        self._orig = {
            "process_sequence": result_db.process_sequence,
            "build_enst_to_refseq_map": result_db.build_enst_to_refseq_map,
            "_load_gene_info": result_db._load_gene_info,
        }
        result_db.process_sequence = _fake_process_sequence
        result_db.build_enst_to_refseq_map = lambda *a, **k: {}
        result_db._load_gene_info = lambda ids: dict((g, LABELS[g]) for g in ids if g in LABELS)

    def tearDown(self):
        for name, fn in self._orig.items():
            setattr(result_db, name, fn)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def build(self, headers=(WT_HEADER, VAR)):
        _write_params(self.output_dir, headers)
        return result_db.ensure_db(self.output_dir)


class SequenceIdTests(unittest.TestCase):
    def test_variant_description_and_label(self):
        self.assertEqual(variant_description(WT), "WT")
        self.assertEqual(variant_description(VAR), "modified 6:U|A")
        self.assertEqual(mirna_node_label(VAR), "hsa-miR-181d-5p (modified 6:U|A)")
        self.assertEqual(mirna_node_label("hsa-miR-21-5p,8:A>U,-4|0,modified_shifted"),
                         "hsa-miR-21-5p (modified+shifted 8:A>U,-4|0)")
        self.assertFalse(is_variant(WT))
        self.assertTrue(is_variant(VAR))

    def test_record_from_parameters_entry(self):
        rec = sequence_record({"header": VAR, "sequence": "ACGU"})
        self.assertEqual(rec["sequence"], VAR)
        self.assertEqual(rec["base_id"], WT)
        self.assertEqual(rec["is_variant"], 1)
        self.assertEqual(rec["seq"], "ACGU")
        wt = sequence_record({"header": WT_HEADER, "sequence": "ACGU"})
        self.assertEqual(wt["sequence"], WT)
        self.assertEqual(wt["variant"], "WT")
        self.assertEqual(wt["is_variant"], 0)

    def test_network_keeps_its_historical_names(self):
        self.assertIs(network_results._mirna_node_id, mirna_node_id)
        self.assertEqual(network_results._mirna_node_label("hsa-miR-21-5p,-7|1,shifted"),
                         "hsa-miR-21-5p (shifted -7|1)")


class BuildTests(_StubbedParsers):
    def test_schema_version_and_sequences_table(self):
        db = self.build()
        self.assertEqual(result_db.schema_version(db), result_db.SCHEMA_VERSION)
        seqs = result_db.list_sequences(db)
        self.assertEqual([s["sequence"] for s in seqs], [WT, VAR])
        self.assertEqual(seqs[0]["is_variant"], False)
        self.assertEqual(seqs[1]["label"], "hsa-miR-181d-5p (modified 6:U|A)")
        self.assertEqual(seqs[1]["header"], VAR)

    def test_hits_are_stored_per_record(self):
        db = self.build()
        self.assertEqual(_tools(db, WT)["NM_021038"], {"miRanda", "PITA", "RNAhybrid"})
        self.assertEqual(_tools(db, VAR)["NM_021038"], {"RNAhybrid"})
        self.assertNotIn("NM_004456", _tools(db, WT))
        self.assertNotIn("NM_002467", _tools(db, VAR))

    def test_skipped_variant_leaves_one_record(self):
        # A modification the runner could not apply never reaches the FASTA,
        # so the job runs the WT only: the selector must show exactly that.
        db = self.build(headers=(WT_HEADER,))
        self.assertEqual([s["sequence"] for s in result_db.list_sequences(db)], [WT])

    def test_version_1_database_is_rebuilt(self):
        _write_params(self.output_dir, (WT_HEADER, VAR))
        db = os.path.join(self.output_dir, result_db.DB_FILENAME)
        conn = sqlite3.connect(db)
        conn.executescript("""
            CREATE TABLE gene_tools (gene_id TEXT, tool TEXT, PRIMARY KEY (gene_id, tool));
            CREATE TABLE gene_info (gene_id TEXT PRIMARY KEY, gene_label TEXT, gene_name TEXT);
            INSERT INTO gene_tools VALUES ('NM_021038', 'miRanda');
        """)
        conn.close()
        self.assertEqual(result_db.schema_version(db), 0)
        result_db.ensure_db(self.output_dir)
        self.assertEqual(result_db.schema_version(db), result_db.SCHEMA_VERSION)
        self.assertEqual(len(result_db.list_sequences(db)), 2)

    def test_current_database_is_not_rebuilt(self):
        db = self.build()
        before = os.stat(db).st_ino
        result_db.ensure_db(self.output_dir)
        self.assertEqual(os.stat(db).st_ino, before)

    def test_predictions_csv(self):
        self.build()
        path = os.path.join(self.output_dir, result_db.PREDICTIONS_CSV)
        with open(path) as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 10)   # 6 WT edges + 4 variant edges
        mbnl1 = [(r["sequence"], r["tool"]) for r in rows if r["gene_id"] == "NM_021038"]
        self.assertEqual(sorted(mbnl1), [(WT, "PITA"), (WT, "RNAhybrid"), (WT, "miRanda"),
                                         (VAR, "RNAhybrid")])
        var_rows = [r for r in rows if r["sequence"] == VAR]
        self.assertEqual(var_rows[0]["mirna"], WT)
        self.assertEqual(var_rows[0]["variant"], "modified 6:U|A")
        self.assertEqual(var_rows[0]["gene_label"], "TP53")


class QueryTests(_StubbedParsers):
    def test_default_view_is_wt_only(self):
        db = self.build()
        res = result_db.query_genes(db, number=100)
        self.assertEqual(res["selected"], [WT])
        self.assertEqual(res["total_genes"], 3)
        by_id = dict((g["gene_id"], g) for g in res["genes"])
        self.assertEqual(by_id["NM_021038"]["tool_count"], 3)
        self.assertNotIn("NM_004456", by_id)
        self.assertEqual(res["venn"]["sets"], {"miRanda": 2, "PITA": 1, "RNAhybrid": 3})

    def test_wt_view_equals_wt_only_job(self):
        # Denise Cilia's check: the WT view of a WT+variant job must equal the
        # result of the same job run without the variant.
        with_variant = result_db.query_genes(self.build(), number=100)
        shutil.rmtree(self.output_dir)
        os.makedirs(self.output_dir)
        wt_only = result_db.query_genes(self.build(headers=(WT_HEADER,)), number=100)
        self.assertEqual(with_variant["genes"], wt_only["genes"])
        self.assertEqual(with_variant["venn"], wt_only["venn"])

    def test_variant_view(self):
        db = self.build()
        res = result_db.query_genes(db, sequences=[VAR], number=100)
        self.assertEqual(res["total_genes"], 3)
        by_id = dict((g["gene_id"], g) for g in res["genes"])
        self.assertEqual(by_id["NM_021038"]["tools"], ["RNAhybrid"])
        self.assertEqual(by_id["NM_004456"]["tools"], ["miRanda"])
        self.assertNotIn("PITA", res["venn"]["sets"])

    def test_union_counts_a_tool_once_per_gene(self):
        db = self.build()
        res = result_db.query_genes(db, sequences=[WT, VAR], number=100)
        self.assertEqual(res["total_genes"], 4)
        by_id = dict((g["gene_id"], g) for g in res["genes"])
        # TP53 is hit by miRanda and RNAhybrid in both records: 2 tools, not 4.
        self.assertEqual(by_id["NM_000546"]["tool_count"], 2)
        self.assertEqual(by_id["NM_000546"]["tools"], ["RNAhybrid", "miRanda"])

    def test_keyword_filter_respects_sequence(self):
        db = self.build()
        res = result_db.query_genes(db, keyword="EZH2", number=100)
        self.assertEqual(res["total"], 0)
        res = result_db.query_genes(db, keyword="EZH2", sequences=[VAR], number=100)
        self.assertEqual(res["total"], 1)
        # Tool keyword: PITA hits only exist for the WT.
        res = result_db.query_genes(db, keyword="PITA", sequences=[VAR], number=100)
        self.assertEqual(res["total"], 0)
        # Set sizes are reported for every tool the record used, all empty here.
        self.assertEqual(res["venn"]["sets"], {"miRanda": 0, "RNAhybrid": 0})
        self.assertEqual(res["venn"]["combinations"], [])

    def test_gene_labels_default_to_wt(self):
        db = self.build()
        self.assertEqual(result_db.gene_labels(db), ["MBNL1", "MYC", "TP53"])
        self.assertEqual(result_db.gene_labels(db, sequences=[VAR]), ["EZH2", "MBNL1", "TP53"])

    def test_version_1_database_still_queryable(self):
        # gene_labels() and query_genes() are also called on databases the
        # tests (and old deployments) built by hand without a sequences table.
        db = os.path.join(self.tmp, "legacy.db")
        conn = sqlite3.connect(db)
        conn.executescript("""
            CREATE TABLE gene_tools (gene_id TEXT, tool TEXT, PRIMARY KEY (gene_id, tool));
            CREATE TABLE gene_info (gene_id TEXT PRIMARY KEY, gene_label TEXT, gene_name TEXT);
            INSERT INTO gene_tools VALUES ('NM_021038', 'miRanda');
            INSERT INTO gene_tools VALUES ('NM_021038', 'PITA');
            INSERT INTO gene_info VALUES ('NM_021038', 'MBNL1', NULL);
        """)
        conn.close()
        res = result_db.query_genes(db)
        self.assertEqual(res["selected"], [])
        self.assertEqual(res["genes"][0]["tool_count"], 2)
        self.assertEqual(result_db.list_sequences(db), [])
        self.assertEqual(result_db.gene_labels(db), ["MBNL1"])


class CompareTests(_StubbedParsers):
    def test_classification(self):
        db = self.build()
        res = result_db.compare_genes(db, WT, VAR, number=100)
        by_id = dict((g["gene_id"], g) for g in res["genes"])
        self.assertEqual(by_id["NM_021038"]["status"], "changed")
        self.assertEqual(by_id["NM_021038"]["lost"], ["PITA", "miRanda"])
        self.assertEqual(by_id["NM_021038"]["gained"], [])
        self.assertEqual(by_id["NM_000546"]["status"], "shared")
        self.assertEqual(by_id["NM_002467"]["status"], "lost")
        self.assertEqual(by_id["NM_002467"]["variant_tools"], [])
        self.assertEqual(by_id["NM_004456"]["status"], "gained")
        self.assertEqual(by_id["NM_004456"]["wt_tools"], [])
        self.assertEqual(res["totals"], {"changed": 1, "lost": 1, "gained": 1, "shared": 1})
        self.assertEqual(res["tools"], {
            "miRanda":   {"gained": 1, "lost": 1},
            "PITA":      {"gained": 0, "lost": 1},
            "RNAhybrid": {"gained": 0, "lost": 1},
        })
        self.assertEqual(res["total_genes"], 4)

    def test_default_order_puts_changes_first(self):
        db = self.build()
        res = result_db.compare_genes(db, WT, VAR, number=100)
        self.assertEqual([g["status"] for g in res["genes"]],
                         ["changed", "lost", "gained", "shared"])
        res = result_db.compare_genes(db, WT, VAR, number=100, order="asc")
        self.assertEqual(res["genes"][0]["status"], "shared")

    def test_status_filter_keyword_and_paging(self):
        db = self.build()
        res = result_db.compare_genes(db, WT, VAR, status="lost")
        self.assertEqual([g["gene_label"] for g in res["genes"]], ["MYC"])
        self.assertEqual(res["total"], 1)
        self.assertEqual(res["total_genes"], 4)          # unfiltered
        self.assertEqual(res["totals"]["changed"], 1)    # unfiltered
        res = result_db.compare_genes(db, WT, VAR, keyword="mbnl")
        self.assertEqual([g["gene_id"] for g in res["genes"]], ["NM_021038"])
        res = result_db.compare_genes(db, WT, VAR, sort_by="gene_label", order="asc",
                                      offset=1, number=2)
        self.assertEqual([g["gene_label"] for g in res["genes"]], ["MBNL1", "MYC"])

    def test_unknown_sequence_compares_as_empty(self):
        db = self.build()
        res = result_db.compare_genes(db, WT, "hsa-miR-181d-5p,nope,shifted", number=100)
        self.assertEqual(res["totals"], {"changed": 0, "lost": 3, "gained": 0, "shared": 0})


class LncrnaBuildTests(unittest.TestCase):
    """The lncRNA builder keeps hits per record too, through the shared writer."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="isotar_lnc_")
        self.output_dir = os.path.join(self.tmp, "output")
        for d in ("miRanda", "PITA"):
            os.makedirs(os.path.join(self.output_dir, d))
        _write_params(self.output_dir, (WT_HEADER, VAR))
        for header in (WT_HEADER, VAR):
            open(os.path.join(self.output_dir, "miRanda",
                              "{}_miRanda_results.txt".format(header)), "w").close()
            open(os.path.join(self.output_dir, "PITA",
                              "{}_PITA_results.tab".format(header)), "w").close()

        def fake(tool):
            def parse(path, result_dict, **_kw):
                header = os.path.basename(path).rsplit("_{}_".format(tool), 1)[0]
                hits = {
                    (WT_HEADER, "miRanda"): ["ENST00000534336"],
                    (WT_HEADER, "PITA"):    ["ENST00000534336", "ENST00000580438"],
                    (VAR, "miRanda"):       ["ENST00000534336"],
                    (VAR, "PITA"):          [],
                }[(header, tool)]
                return {"prediction": {tool: hits}}
            return parse

        self._orig = {
            "parseMirandaResults": lncrna_results.parseMirandaResults,
            "parsePITAResults": lncrna_results.parsePITAResults,
        }
        lncrna_results.parseMirandaResults = fake("miRanda")
        lncrna_results.parsePITAResults = fake("PITA")

    def tearDown(self):
        for name, fn in self._orig.items():
            setattr(lncrna_results, name, fn)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_per_record_rows_and_transcript_labels(self):
        db = lncrna_results.ensure_lncrna_db(self.output_dir)
        self.assertEqual(result_db.schema_version(db), result_db.SCHEMA_VERSION)
        self.assertEqual(_tools(db, WT)["ENST00000534336"], {"miRanda", "PITA"})
        self.assertEqual(_tools(db, VAR)["ENST00000534336"], {"miRanda"})
        self.assertNotIn("ENST00000580438", _tools(db, VAR))
        res = result_db.query_genes(db, number=100)
        self.assertEqual(res["selected"], [WT])
        self.assertEqual(res["total_genes"], 2)
        self.assertEqual(res["genes"][0]["gene_label"], res["genes"][0]["gene_id"])
        cmp_ = result_db.compare_genes(db, WT, VAR, number=100)
        self.assertEqual(cmp_["totals"], {"changed": 1, "lost": 1, "gained": 0, "shared": 0})
        self.assertTrue(os.path.exists(
            os.path.join(self.output_dir, lncrna_results.LNCRNA_PREDICTIONS_CSV)))

    def test_network_pool_reads_the_same_edges(self):
        edges = network_results._collect_lncrna_edges(self.output_dir)
        self.assertEqual(edges[WT]["ENST00000534336"], {"miRanda", "PITA"})
        self.assertEqual(edges[VAR]["ENST00000534336"], {"miRanda"})
        self.assertEqual(network_results._collect_lncrna_edges(os.path.join(self.tmp, "nope")), {})


if __name__ == "__main__":
    unittest.main()
