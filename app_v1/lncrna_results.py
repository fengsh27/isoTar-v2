"""Lightweight result reader for the miR-LncRNA workflow.

The gene flow (result_db._build_db) aggregates predictions by protein-coding
gene: it converts target IDs to RefSeq and looks up HGNC symbols/names from
reference_mapping.db. None of that applies to lncRNA targets, whose hits are
Ensembl lncRNA transcripts (ENST...) with no RefSeq/gene-symbol identity.

This module builds the same SQLite schema (sequences/gene_tools/gene_info, see
result_db.write_db) so the existing result_db.query_genes() pagination + Venn
logic can be reused unchanged -- but it aggregates by lncRNA transcript ID and
performs NO gene mapping:

  * only the five lncRNA-compatible tools are read (miRanda, RNAhybrid, miRmap,
    DMISO, PITA); TargetScan is rejected upstream for this workflow (it ignores
    the target FASTA and reads its own precomputed 3' UTR datasets).
  * each tool's raw output file is read with parse_result's per-tool parser,
    which extracts the ENST transcript and applies the same biophysical hit
    filters (seed match, dG/ddG cutoffs) as the gene flow.
  * gene_info.gene_label is set to the transcript ID and gene_name to NULL --
    no reference_mapping.db lookup.

Hits are kept per miRNA record (WT and each variant separately), exactly as the
gene flow does; the network builder reuses collect_lncrna_edges() for its
lncRNA pool.
"""

import os

from app_v1.parse_result import (
    read_sequences_from_json,
    parseMirandaResults,
    parseRnahybridResults,
    parseMirmapResults,
    parseDMISOResults,
    parsePITAResults,
    _extract_lncrna_transcript_id,
)
from app_v1.result_db import SCHEMA_VERSION, schema_version, write_db
from app_v1.sequence_ids import sequence_record

LNCRNA_DB_FILENAME = "lncrna_result.db"
LNCRNA_PREDICTIONS_CSV = "lncrna_predictions.csv"


def collect_lncrna_edges(output_dir):
    """Parse every record's raw lncRNA tool output in output_dir.

    Returns (sequence_rows, edges):
        sequence_rows: one sequence_record() per record in
                       mirna_prediction_parameters.json, in file order
        edges:         sequence id -> {transcript_id -> set(tool names)}
    """
    json_file = os.path.join(output_dir, "mirna_prediction_parameters.json")
    sequences = read_sequences_from_json(json_file)
    ext = _extract_lncrna_transcript_id

    sequence_rows = []
    edges = {}
    for sequence in sequences:
        header = sequence["header"]
        row = sequence_record(sequence)
        sequence_rows.append(row)
        targets_of = edges.setdefault(row["sequence"], {})

        def _add(transcript_id, tool):
            if transcript_id:
                targets_of.setdefault(transcript_id, set()).add(tool)

        p_miranda = os.path.join(output_dir, "miRanda", "{}_miRanda_results.txt".format(header))
        if os.path.exists(p_miranda):
            r = parseMirandaResults(p_miranda, {}, id_extractor=ext)
            for tid in r["prediction"].get("miRanda", []):
                _add(tid, "miRanda")

        p_rnahybrid = os.path.join(output_dir, "RNAhybrid", "{}_RNAhybrid_results.txt".format(header))
        if os.path.exists(p_rnahybrid):
            r = parseRnahybridResults(p_rnahybrid, {}, id_extractor=ext)
            for tid in r["prediction"].get("RNAhybrid", []):
                _add(tid, "RNAhybrid")

        p_mirmap = os.path.join(output_dir, "miRmap", "{}_miRmap_results.txt".format(header))
        if os.path.exists(p_mirmap):
            r = parseMirmapResults(p_mirmap, {}, id_extractor=ext)
            for tid in r["prediction"].get("miRmap", []):
                _add(tid, "miRmap")

        p_dmiso = os.path.join(output_dir, "DMISO", "{}_DMISO_results.txt".format(header))
        if os.path.exists(p_dmiso):
            r = parseDMISOResults(
                p_dmiso, {},
                mirna_sequence=sequence.get("sequence"),
                id_extractor=ext,
            )
            for tid in r["prediction"].get("DMISO", []):
                _add(tid, "DMISO")

        p_pita = os.path.join(output_dir, "PITA", "{}_PITA_results.tab".format(header))
        if os.path.exists(p_pita):
            r = parsePITAResults(p_pita, {}, id_extractor=ext)
            for tid in r["prediction"].get("PITA", []):
                _add(tid, "PITA")

    return sequence_rows, edges


def _build_lncrna_db(output_dir, db_path):
    """Parse lncRNA predictions and populate a fresh SQLite database whose schema
    matches result_db so query_genes() can read it. gene_id == gene_label ==
    transcript ID; gene_name is NULL (no gene mapping for lncRNA).

    Writes to a temp file then renames atomically (partial writes never visible).
    """
    sequence_rows, edges = collect_lncrna_edges(output_dir)

    transcript_ids = set()
    for targets_of in edges.values():
        transcript_ids.update(targets_of.keys())
    # gene_label = transcript ID so the existing table/sort/filter render the
    # lncRNA accession; gene_name stays NULL (no symbol identity).
    gene_info = dict((tid, (tid, None)) for tid in transcript_ids)

    write_db(output_dir, db_path, sequence_rows, edges, gene_info,
             csv_name=LNCRNA_PREDICTIONS_CSV)


def ensure_lncrna_db(output_dir):
    """Return path to lncrna_result.db, building it first if it does not exist
    or was written by an older schema (see result_db.SCHEMA_VERSION)."""
    db_path = os.path.join(output_dir, LNCRNA_DB_FILENAME)
    if not os.path.exists(db_path) or schema_version(db_path) < SCHEMA_VERSION:
        _build_lncrna_db(output_dir, db_path)
    return db_path
