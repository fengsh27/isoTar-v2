"""Builder for the miR-Network workflow (tripartite gene <-> miRNA <-> lncRNA).

A mir-network job runs BOTH prediction pools over the same miRNA list:

    output/gene/    -- 3' UTR (protein-coding) pool, all selected tools
    output/lncrna/  -- lncRNA transcript pool, TargetScan dropped

Unlike result_db / lncrna_results -- which aggregate by target and throw away
which miRNA produced each hit -- the network needs the (miRNA -> target) edges.
Those are recoverable because mirna_predicting.py writes one output file per
(miRNA, tool) named by the miRNA header (e.g. "hsa-miR-495-3p,WT_miRanda_..."),
so parsing per sequence keeps the miRNA attribution.

The graph is assembled in two layers:

  * file parsing  -- _collect_pool_edges() reuses parse_result's per-tool parsers
    (with the gene flow's ENST->RefSeq mapping on the gene side and the lenient
    lncRNA transcript extractor on the lncRNA side) to produce, per pool, a map
    mirna_id -> {target_id -> set(tools)}.

  * graph shaping -- compute_network() is a PURE function over those two maps
    (no filesystem) that applies the ceRNA "co-target bridge" semantics: for a
    (gene, lncRNA) pair, the bridge miRNAs are those predicted to target BOTH
    members. With user pairs it filters to them; without pairs it falls back to
    a bounded discovery mode over the top-degree genes/lncRNAs.

Parsed edges are cached in network.db (same lazy-build pattern as result.db) so
re-querying with different thresholds/pairs does not re-parse the tool output.
"""

import json
import os
import re
import sqlite3
import tempfile

from app_v1.parse_result import (
    read_sequences_from_json,
    process_sequence,
    build_enst_to_refseq_map,
    parseMirandaResults,
    parseRnahybridResults,
    parseMirmapResults,
    parseDMISOResults,
    parsePITAResults,
    _extract_lncrna_transcript_id,
)
from app_v1.result_db import _load_gene_info

NETWORK_DB_FILENAME = "network.db"

# Discovery-mode caps (no user pairs): bound the combinatorial gene x lncRNA
# bridge search to the highest-degree targets so the graph stays renderable.
DEFAULT_TOP_GENES = 50
DEFAULT_TOP_LNCRNA = 50
MAX_DISCOVERY_PAIRS = 500

_VERSION_SUFFIX_RE = re.compile(r"\.\d+$")


def strip_lncrna_version(transcript_id):
    """Strip a trailing Ensembl/FlyBase '.<version>' from a lncRNA transcript ID.

    Mirrors _extract_lncrna_transcript_id's version handling so user-supplied
    pair lncRNAs match the parsed hit IDs. WormBase isoform names (Y51H4A.27)
    are left intact -- only a purely numeric Ensembl-style suffix is removed.
    """
    if not transcript_id:
        return transcript_id
    t = transcript_id.strip()
    extracted = _extract_lncrna_transcript_id(t)
    return extracted or t


# ---------------------------------------------------------------------------
# File parsing -> per-pool (mirna -> target -> tools) maps
# ---------------------------------------------------------------------------

def _mirna_node_id(header):
    """Node id for a miRNA record, from the runner header "<id>,<variant...>".

    A WT record keeps the plain base id (so a miRNA with no shift/modification
    is a single node, exactly as before). A modified/shifted variant keeps its
    full header, so it becomes a node distinct from its WT baseline -- letting
    the graph show how the modification redirects targeting.
        "hsa-miR-21-5p,WT"            -> "hsa-miR-21-5p"
        "hsa-miR-21-5p,-7|1,shifted"  -> "hsa-miR-21-5p,-7|1,shifted"
    """
    parts = header.split(",")
    if len(parts) <= 1 or parts[-1] == "WT":
        return parts[0]
    return header


def _mirna_base_id(node_id):
    """Base miRNA id a (possibly variant) node belongs to, for grouping."""
    return node_id.split(",")[0]


def _mirna_node_label(node_id):
    """Human-readable label, e.g. 'hsa-miR-21-5p (shifted -7|1)'. WT/plain nodes
    keep the bare id."""
    parts = node_id.split(",")
    if len(parts) <= 1:
        return parts[0]
    base, vtype = parts[0], parts[-1].replace("_", "+")
    detail = ",".join(parts[1:-1])
    return "{} ({} {})".format(base, vtype, detail).strip()


def _collect_gene_edges(pool_dir, enst_to_refseq):
    """Parse the gene (3' UTR) pool into mirna_id -> {refseq -> set(tools)}.

    Reuses parse_result.process_sequence (the gene flow's parser), which yields
    RefSeq IDs per tool with TargetScan ENSTs mapped to RefSeq. No target filter
    is applied -- the ceRNA pairs do the filtering downstream."""
    edges = {}
    json_file = os.path.join(pool_dir, "mirna_prediction_parameters.json")
    if not os.path.exists(json_file):
        return edges
    for sequence in read_sequences_from_json(json_file):
        mirna_id = _mirna_node_id(sequence["header"])
        targets = edges.setdefault(mirna_id, {})
        results = process_sequence(
            sequence, pool_dir, enst_to_refseq=enst_to_refseq, targets=None,
        )
        for tool, gene_ids in results.get("prediction", {}).items():
            for gene_id in gene_ids:
                if gene_id:
                    targets.setdefault(gene_id, set()).add(tool)
    return edges


def _collect_lncrna_edges(pool_dir):
    """Parse the lncRNA pool into mirna_id -> {transcript_id -> set(tools)}.

    Uses the five lncRNA-compatible tools with the lenient lncRNA transcript
    extractor (no gene mapping), mirroring lncrna_results but keeping the miRNA."""
    edges = {}
    json_file = os.path.join(pool_dir, "mirna_prediction_parameters.json")
    if not os.path.exists(json_file):
        return edges
    ext = _extract_lncrna_transcript_id
    for sequence in read_sequences_from_json(json_file):
        header = sequence["header"]
        mirna_id = _mirna_node_id(header)
        targets = edges.setdefault(mirna_id, {})

        def _add(tid, tool):
            if tid:
                targets.setdefault(tid, set()).add(tool)

        p = os.path.join(pool_dir, "miRanda", "{}_miRanda_results.txt".format(header))
        if os.path.exists(p):
            r = parseMirandaResults(p, {}, id_extractor=ext)
            for tid in r["prediction"].get("miRanda", []):
                _add(tid, "miRanda")

        p = os.path.join(pool_dir, "RNAhybrid", "{}_RNAhybrid_results.txt".format(header))
        if os.path.exists(p):
            r = parseRnahybridResults(p, {}, id_extractor=ext)
            for tid in r["prediction"].get("RNAhybrid", []):
                _add(tid, "RNAhybrid")

        p = os.path.join(pool_dir, "miRmap", "{}_miRmap_results.txt".format(header))
        if os.path.exists(p):
            r = parseMirmapResults(p, {}, id_extractor=ext)
            for tid in r["prediction"].get("miRmap", []):
                _add(tid, "miRmap")

        p = os.path.join(pool_dir, "DMISO", "{}_DMISO_results.txt".format(header))
        if os.path.exists(p):
            r = parseDMISOResults(
                p, {}, mirna_sequence=sequence.get("sequence"), id_extractor=ext,
            )
            for tid in r["prediction"].get("DMISO", []):
                _add(tid, "DMISO")

        p = os.path.join(pool_dir, "PITA", "{}_PITA_results.tab".format(header))
        if os.path.exists(p):
            r = parsePITAResults(p, {}, id_extractor=ext)
            for tid in r["prediction"].get("PITA", []):
                _add(tid, "PITA")

    return edges


# ---------------------------------------------------------------------------
# network.db cache (parsed edges)
# ---------------------------------------------------------------------------

def _build_network_db(output_dir, db_path):
    """Parse both pools and persist the (mirna, target, target_type, tool) edge
    set plus per-target labels. Atomic temp-file-then-rename like result_db."""
    gene_dir = os.path.join(output_dir, "gene")
    lncrna_dir = os.path.join(output_dir, "lncrna")

    enst_to_refseq = build_enst_to_refseq_map()
    gene_edges = _collect_gene_edges(gene_dir, enst_to_refseq)
    lncrna_edges = _collect_lncrna_edges(lncrna_dir)

    gene_ids = set()
    for targets in gene_edges.values():
        gene_ids.update(targets.keys())
    gene_info = _load_gene_info(gene_ids)

    tmp_fd, tmp_path = tempfile.mkstemp(dir=output_dir, suffix=".db.tmp")
    os.close(tmp_fd)
    try:
        conn = sqlite3.connect(tmp_path)
        try:
            c = conn.cursor()
            c.execute("""
                CREATE TABLE edges (
                    mirna       TEXT NOT NULL,
                    target      TEXT NOT NULL,
                    target_type TEXT NOT NULL,
                    tool        TEXT NOT NULL,
                    PRIMARY KEY (mirna, target, target_type, tool)
                )
            """)
            c.execute("""
                CREATE TABLE labels (
                    target      TEXT NOT NULL,
                    target_type TEXT NOT NULL,
                    label       TEXT,
                    name        TEXT,
                    PRIMARY KEY (target, target_type)
                )
            """)
            c.executemany(
                "INSERT OR IGNORE INTO edges (mirna, target, target_type, tool) VALUES (?, ?, 'gene', ?)",
                ((mirna, target, tool)
                 for mirna, targets in gene_edges.items()
                 for target, tools in targets.items()
                 for tool in tools),
            )
            c.executemany(
                "INSERT OR IGNORE INTO edges (mirna, target, target_type, tool) VALUES (?, ?, 'lncrna', ?)",
                ((mirna, target, tool)
                 for mirna, targets in lncrna_edges.items()
                 for target, tools in targets.items()
                 for tool in tools),
            )
            c.executemany(
                "INSERT OR IGNORE INTO labels (target, target_type, label, name) VALUES (?, 'gene', ?, ?)",
                ((gid, gene_info.get(gid, (None, None))[0] or gid,
                  gene_info.get(gid, (None, None))[1]) for gid in gene_ids),
            )
            conn.commit()
        finally:
            conn.close()
        os.rename(tmp_path, db_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


def ensure_network_db(output_dir):
    """Return path to network.db, building it from both pools if it is absent."""
    db_path = os.path.join(output_dir, NETWORK_DB_FILENAME)
    if not os.path.exists(db_path):
        _build_network_db(output_dir, db_path)
    return db_path


def _load_edges(db_path):
    """Read network.db into (gene_edges, lncrna_edges, labels):
        gene_edges/lncrna_edges: target -> {mirna -> set(tools)}
        labels:                  (target, target_type) -> (label, name)
    Keyed by target (not miRNA) because the bridge search intersects the miRNA
    sets of a gene and a lncRNA."""
    gene_edges = {}
    lncrna_edges = {}
    labels = {}
    conn = sqlite3.connect(db_path)
    try:
        c = conn.cursor()
        c.execute("SELECT mirna, target, target_type, tool FROM edges")
        for mirna, target, target_type, tool in c.fetchall():
            bucket = gene_edges if target_type == "gene" else lncrna_edges
            bucket.setdefault(target, {}).setdefault(mirna, set()).add(tool)
        c.execute("SELECT target, target_type, label, name FROM labels")
        for target, target_type, label, name in c.fetchall():
            labels[(target, target_type)] = (label, name)
    finally:
        conn.close()
    return gene_edges, lncrna_edges, labels


# ---------------------------------------------------------------------------
# Pure graph shaping (no filesystem) -- unit-testable in isolation
# ---------------------------------------------------------------------------

def compute_network(gene_edges, lncrna_edges, pairs=None, labels=None,
                    top_genes=DEFAULT_TOP_GENES, top_lncrna=DEFAULT_TOP_LNCRNA,
                    max_pairs=MAX_DISCOVERY_PAIRS):
    """Build the tripartite network JSON from per-target edge maps.

    Args:
        gene_edges:   gene_id -> {mirna -> set(tools)}
        lncrna_edges: lncrna_id -> {mirna -> set(tools)}
        pairs:        optional list of {"gene": <refseq or symbol input>,
                      "gene_refseqs": [refseq, ...], "lncrna": <transcript>}.
                      When present -> ceRNA filter mode (only these pairs).
                      When empty/None -> discovery mode over top-degree targets.
        labels:       (target, target_type) -> (label, name) for display.
        top_genes/top_lncrna/max_pairs: discovery-mode bounds.

    Returns a dict {mode, nodes, edges, pairs, summary}. Every gene-miRNA and
    miRNA-lncRNA edge carries its full tool list so the client can threshold by
    consensus k without a server round-trip.
    """
    labels = labels or {}

    if pairs:
        mode = "pairs"
        candidate_pairs = []
        for p in pairs:
            lncrna = strip_lncrna_version(p.get("lncrna"))
            refseqs = p.get("gene_refseqs") or ([p["gene"]] if p.get("gene") else [])
            for gene in refseqs:
                candidate_pairs.append((gene, lncrna, p.get("gene")))
        truncated = False
    else:
        mode = "discovery"
        # Rank by number of distinct miRNAs targeting the node (bridge potential).
        top_g = sorted(gene_edges, key=lambda g: len(gene_edges[g]), reverse=True)[:top_genes]
        top_l = sorted(lncrna_edges, key=lambda l: len(lncrna_edges[l]), reverse=True)[:top_lncrna]
        candidate_pairs = [(g, l, None) for g in top_g for l in top_l]
        truncated = len(candidate_pairs) > max_pairs
        candidate_pairs = candidate_pairs[:max_pairs]

    gene_nodes = {}    # gene_id -> node
    lncrna_nodes = {}  # lncrna_id -> node
    mirna_nodes = {}   # mirna_id -> node
    gene_mirna = {}    # (gene_id, mirna) -> set(tools)
    mirna_lncrna = {}  # (mirna, lncrna_id) -> set(tools)
    kept_pairs = []

    for gene, lncrna, gene_input in candidate_pairs:
        g_mirnas = gene_edges.get(gene)
        l_mirnas = lncrna_edges.get(lncrna)
        if not g_mirnas or not l_mirnas:
            continue
        bridges = set(g_mirnas) & set(l_mirnas)
        if not bridges:
            continue

        g_label = labels.get((gene, "gene"), (None, None))
        gene_nodes.setdefault(gene, {
            "id": gene, "type": "gene",
            "label": g_label[0] or gene, "name": g_label[1],
        })
        lncrna_nodes.setdefault(lncrna, {
            "id": lncrna, "type": "lncrna", "label": lncrna, "name": None,
        })
        for m in bridges:
            mirna_nodes.setdefault(m, {
                "id": m, "type": "mirna",
                "label": _mirna_node_label(m), "base": _mirna_base_id(m),
            })
            gene_mirna.setdefault((gene, m), set()).update(g_mirnas[m])
            mirna_lncrna.setdefault((m, lncrna), set()).update(l_mirnas[m])

        kept_pairs.append({
            "gene": gene,
            "gene_input": gene_input,
            "gene_label": g_label[0] or gene,
            "lncrna": lncrna,
            "bridge_mirnas": sorted(bridges),
        })

    edges = []
    for (gene, m), tools in gene_mirna.items():
        edges.append({
            "source": gene, "target": m, "side": "gene",
            "tools": sorted(tools), "tool_count": len(tools),
        })
    for (m, lncrna), tools in mirna_lncrna.items():
        edges.append({
            "source": m, "target": lncrna, "side": "lncrna",
            "tools": sorted(tools), "tool_count": len(tools),
        })

    nodes = list(gene_nodes.values()) + list(mirna_nodes.values()) + list(lncrna_nodes.values())
    summary = {
        "mode": mode,
        "gene_count": len(gene_nodes),
        "mirna_count": len(mirna_nodes),
        "lncrna_count": len(lncrna_nodes),
        "edge_count": len(edges),
        "pair_count": len(kept_pairs),
        "truncated": truncated,
    }
    return {"mode": mode, "nodes": nodes, "edges": edges,
            "pairs": kept_pairs, "summary": summary}


def build_network(output_dir, pairs=None, top_genes=DEFAULT_TOP_GENES,
                  top_lncrna=DEFAULT_TOP_LNCRNA):
    """Build (and cache) network.db, then shape the tripartite graph JSON."""
    db_path = ensure_network_db(output_dir)
    gene_edges, lncrna_edges, labels = _load_edges(db_path)
    return compute_network(
        gene_edges, lncrna_edges, pairs=pairs, labels=labels,
        top_genes=top_genes, top_lncrna=top_lncrna,
    )


def load_pairs(job_dir):
    """Read resolved pairs.json (written at submit time) from the job dir.
    Returns the list of pair dicts, or [] if absent/unreadable."""
    path = os.path.join(job_dir, "pairs.json")
    if not os.path.exists(path):
        return []
    try:
        with open(path, "r") as f:
            data = json.load(f)
    except (ValueError, IOError):
        return []
    return data.get("pairs", []) if isinstance(data, dict) else []
