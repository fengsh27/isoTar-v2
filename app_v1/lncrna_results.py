"""Lightweight result reader for the miR-LncRNA workflow.

The gene flow (result_db._build_db) aggregates predictions by protein-coding
gene: it converts target IDs to RefSeq and looks up HGNC symbols/names from
reference_mapping.db. None of that applies to lncRNA targets, whose hits are
Ensembl lncRNA transcripts (ENST...) with no RefSeq/gene-symbol identity.

This module builds the same SQLite schema (gene_tools/gene_info) so the existing
result_db.query_genes() pagination + Venn logic can be reused unchanged -- but it
aggregates by lncRNA transcript ID and performs NO gene mapping:

  * only the four lncRNA-compatible tools are read (miRanda, RNAhybrid, miRmap,
    DMISO); TargetScan/PITA are rejected upstream for this workflow.
  * each tool's raw output file is read with parse_result's per-tool parser,
    which extracts the ENST transcript and applies the same biophysical hit
    filters (seed match, dG/ddG cutoffs) as the gene flow.
  * gene_info.gene_label is set to the transcript ID and gene_name to NULL --
    no reference_mapping.db lookup.
"""

import os
import sqlite3
import tempfile

from app_v1.parse_result import (
    read_sequences_from_json,
    parseMirandaResults,
    parseRnahybridResults,
    parseMirmapResults,
    parseDMISOResults,
    _extract_lncrna_transcript_id,
)

LNCRNA_DB_FILENAME = "lncrna_result.db"


def _collect_transcript_tools(output_dir):
    """Walk the four lncRNA-compatible tool output dirs and return a dict
    transcript_id -> set(tool names) of predicted hits."""
    json_file = os.path.join(output_dir, "mirna_prediction_parameters.json")
    sequences = read_sequences_from_json(json_file)

    transcript_tools = {}

    def _add(transcript_id, tool):
        if not transcript_id:
            return
        transcript_tools.setdefault(transcript_id, set()).add(tool)

    for sequence in sequences:
        header = sequence["header"]

        p_miranda = os.path.join(output_dir, "miRanda", "{}_miRanda_results.txt".format(header))
        if os.path.exists(p_miranda):
            r = parseMirandaResults(p_miranda, {}, id_extractor=_extract_lncrna_transcript_id)
            for tid in r["prediction"].get("miRanda", []):
                _add(tid, "miRanda")

        p_rnahybrid = os.path.join(output_dir, "RNAhybrid", "{}_RNAhybrid_results.txt".format(header))
        if os.path.exists(p_rnahybrid):
            r = parseRnahybridResults(p_rnahybrid, {}, id_extractor=_extract_lncrna_transcript_id)
            for tid in r["prediction"].get("RNAhybrid", []):
                _add(tid, "RNAhybrid")

        p_mirmap = os.path.join(output_dir, "miRmap", "{}_miRmap_results.txt".format(header))
        if os.path.exists(p_mirmap):
            r = parseMirmapResults(p_mirmap, {}, id_extractor=_extract_lncrna_transcript_id)
            for tid in r["prediction"].get("miRmap", []):
                _add(tid, "miRmap")

        p_dmiso = os.path.join(output_dir, "DMISO", "{}_DMISO_results.txt".format(header))
        if os.path.exists(p_dmiso):
            r = parseDMISOResults(
                p_dmiso, {},
                mirna_sequence=sequence.get("sequence"),
                id_extractor=_extract_lncrna_transcript_id,
            )
            for tid in r["prediction"].get("DMISO", []):
                _add(tid, "DMISO")

    return transcript_tools


def _build_lncrna_db(output_dir, db_path):
    """Parse lncRNA predictions and populate a fresh SQLite database whose schema
    matches result_db so query_genes() can read it. gene_id == gene_label ==
    transcript ID; gene_name is NULL (no gene mapping for lncRNA).

    Writes to a temp file then renames atomically (partial writes never visible).
    """
    transcript_tools = _collect_transcript_tools(output_dir)

    tmp_fd, tmp_path = tempfile.mkstemp(dir=output_dir, suffix=".db.tmp")
    os.close(tmp_fd)
    try:
        conn = sqlite3.connect(tmp_path)
        try:
            c = conn.cursor()
            c.execute("""
                CREATE TABLE gene_tools (
                    gene_id TEXT NOT NULL,
                    tool    TEXT NOT NULL,
                    PRIMARY KEY (gene_id, tool)
                )
            """)
            c.execute("""
                CREATE TABLE gene_info (
                    gene_id    TEXT PRIMARY KEY,
                    gene_label TEXT,
                    gene_name  TEXT
                )
            """)
            c.executemany(
                "INSERT OR IGNORE INTO gene_tools (gene_id, tool) VALUES (?, ?)",
                ((transcript_id, tool)
                 for transcript_id, tools in transcript_tools.items()
                 for tool in tools),
            )
            # gene_label = transcript ID so the existing table/sort/filter render
            # the lncRNA accession; gene_name stays NULL (no symbol identity).
            c.executemany(
                "INSERT OR IGNORE INTO gene_info (gene_id, gene_label, gene_name) VALUES (?, ?, NULL)",
                ((transcript_id, transcript_id) for transcript_id in transcript_tools.keys()),
            )
            conn.commit()
        finally:
            conn.close()
        os.rename(tmp_path, db_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


def ensure_lncrna_db(output_dir):
    """Return path to lncrna_result.db, building it first if it does not exist."""
    db_path = os.path.join(output_dir, LNCRNA_DB_FILENAME)
    if not os.path.exists(db_path):
        _build_lncrna_db(output_dir, db_path)
    return db_path
