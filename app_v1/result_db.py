import os
import sqlite3
import sys
import tempfile
from itertools import combinations

# Allow importing parse_result from the v2 package (works in dev and in Docker)
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from app_v1.parse_result import read_sequences_from_json, process_sequence, _extract_transcript_id

DB_FILENAME = "result.db"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_dmiso_file(dmiso_path):
    """Return deduplicated gene IDs from a pre-filtered DMISO results file.

    File format (tab-separated, score already filtered > 0.99 by mirna_predicting.py):
        Target ID\\tTarget Sequence\\tPrediction Score
        hg19_refGene_NM_144722 range=...\\t<seq>\\t0.991...
    """
    gene_ids = set()
    with open(dmiso_path, "r") as f:
        next(f)  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            target_id = line.split("\t")[0]   # first column
            gene_id = _extract_transcript_id(target_id)
            if gene_id:
                gene_ids.add(gene_id)
    return gene_ids


def _build_db(output_dir, db_path):
    """Parse all prediction results and populate a fresh SQLite database.

    Schema:
        gene_tools(gene_id TEXT, tool TEXT, PRIMARY KEY (gene_id, tool))

    Non-DMISO tools are parsed via parse_result.process_sequence().
    DMISO is read directly from output/DMISO/<header>_DMISO_results.txt
    because parseDMISOResults() in parse_result.py is a no-op.

    Writes to a temp file first then renames atomically so a partial write
    is never visible at db_path.
    """
    json_file = os.path.join(output_dir, "mirna_prediction_parameters.json")
    sequences = read_sequences_from_json(json_file)

    # Collect gene → tools mapping across all sequences
    gene_tools = {}   # gene_id -> set of tool names

    def _add(gene_id, tool):
        if gene_id is None:
            return
        if gene_id not in gene_tools:
            gene_tools[gene_id] = set()
        gene_tools[gene_id].add(tool)

    for sequence in sequences:
        # Non-DMISO tools via parse_result
        results = process_sequence(sequence, output_dir)
        for tool, gene_ids in results.get("prediction", {}).items():
            for gene_id in gene_ids:
                _add(gene_id, tool)

        # DMISO — read directly (parseDMISOResults is a no-op)
        dmiso_path = os.path.join(
            output_dir, "DMISO",
            "{}_DMISO_results.txt".format(sequence["header"]),
        )
        if os.path.exists(dmiso_path):
            for gene_id in _parse_dmiso_file(dmiso_path):
                _add(gene_id, "DMISO")

    # Write to a temp file, then rename atomically
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
            c.executemany(
                "INSERT OR IGNORE INTO gene_tools (gene_id, tool) VALUES (?, ?)",
                ((gene_id, tool)
                 for gene_id, tools in gene_tools.items()
                 for tool in tools),
            )
            conn.commit()
        finally:
            conn.close()
        os.rename(tmp_path, db_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


def _venn_stats(cursor, filter_pattern=None):
    """Compute per-tool set sizes and all pairwise/higher-order intersections.

    Intersection semantics: |A∩B| counts genes found by BOTH A and B,
    regardless of other tools — the standard input for Venn diagram libraries.

    filter_pattern: SQL LIKE pattern (e.g. '%ENST%'), or None for no filter.
    """
    cursor.execute("SELECT DISTINCT tool FROM gene_tools ORDER BY tool")
    tools = [row[0] for row in cursor.fetchall()]

    gene_filter_sql = "AND gene_id LIKE ?" if filter_pattern else ""
    gene_filter_arg = [filter_pattern] if filter_pattern else []

    sets = {}
    for tool in tools:
        cursor.execute(
            "SELECT COUNT(DISTINCT gene_id) FROM gene_tools WHERE tool = ? {}".format(
                gene_filter_sql
            ),
            [tool] + gene_filter_arg,
        )
        sets[tool] = cursor.fetchone()[0]

    intersections = {}
    for r in range(2, len(tools) + 1):
        for combo in combinations(tools, r):
            placeholders = ",".join("?" * len(combo))
            cursor.execute(
                """
                SELECT COUNT(*) FROM (
                    SELECT gene_id FROM gene_tools
                    WHERE  tool IN ({}) {}
                    GROUP  BY gene_id
                    HAVING COUNT(DISTINCT tool) = ?
                )
                """.format(placeholders, gene_filter_sql),
                list(combo) + gene_filter_arg + [len(combo)],
            )
            count = cursor.fetchone()[0]
            if count > 0:
                intersections["&".join(combo)] = count

    return {"sets": sets, "intersections": intersections}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def ensure_db(output_dir):
    """Return path to result.db, building it first if it does not exist."""
    db_path = os.path.join(output_dir, DB_FILENAME)
    if not os.path.exists(db_path):
        _build_db(output_dir, db_path)
    return db_path


def query_genes(db_path, sort_by="tool_count", order="desc", offset=0, number=20,
                gene_label=None):
    """Query the gene-tool database and return paginated results + Venn stats.

    Args:
        sort_by:    "tool_count" | "gene_label"
        order:      "asc" | "desc"
        offset:     row offset for pagination
        number:     page size
        gene_label: optional substring to filter gene_id (case-insensitive)

    Returns a dict with keys: total, genes, venn.
    """
    sort_col = "tool_count" if sort_by == "tool_count" else "gene_id"
    sort_dir = "ASC" if order.lower() in ("asc", "ascend") else "DESC"

    filter_pattern = "%{}%".format(gene_label) if gene_label else None
    where_sql  = "WHERE gene_id LIKE ?" if filter_pattern else ""
    where_args = [filter_pattern] if filter_pattern else []

    conn = sqlite3.connect(db_path)
    try:
        c = conn.cursor()

        # Total distinct genes — unfiltered (always the full result set)
        c.execute("SELECT COUNT(DISTINCT gene_id) FROM gene_tools")
        total_genes = c.fetchone()[0]

        # Total distinct genes — respects geneLabel filter
        c.execute(
            "SELECT COUNT(DISTINCT gene_id) FROM gene_tools {}".format(where_sql),
            where_args,
        )
        total = c.fetchone()[0]

        # Paginated gene rows (respects filter)
        c.execute(
            """
            SELECT   gene_id,
                     COUNT(*)               AS tool_count,
                     GROUP_CONCAT(tool, ',') AS tools
            FROM     gene_tools
            {where}
            GROUP BY gene_id
            ORDER BY {sort_col} {sort_dir}, gene_id ASC
            LIMIT  ? OFFSET ?
            """.format(where=where_sql, sort_col=sort_col, sort_dir=sort_dir),
            where_args + [number, offset],
        )
        genes = [
            {
                "gene_id":    row[0],
                "tool_count": row[1],
                "tools":      sorted(row[2].split(",")),
            }
            for row in c.fetchall()
        ]

        venn = _venn_stats(c, filter_pattern=filter_pattern)

    finally:
        conn.close()

    return {"total_genes": total_genes, "total": total, "genes": genes, "venn": venn}


if __name__ == "__main__":
    ensure_db("out")
    res = query_genes("out")
    print(res)
