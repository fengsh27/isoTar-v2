import json
import os
import sqlite3
import sys
import tempfile
from itertools import combinations

# Allow importing parse_result from the v2 package (works in dev and in Docker)
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from app_v1.parse_result import (
    read_sequences_from_json,
    process_sequence,
    build_enst_to_refseq_map,
    load_targets_file,
)

DB_FILENAME = "result.db"

# Path to the reference mapping database (gene_id -> gene_label/gene_name)
REFERENCE_MAPPING_DB = os.environ.get(
    "ISOTAR_REFERENCE_MAPPING_DB",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference_mapping.db"),
)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_gene_info(gene_ids):
    """Look up gene_label and gene_name for a collection of gene IDs.

    Returns a dict: gene_id -> (gene_label, gene_name).
    Missing entries map to (None, None).
    """
    info = {}
    if not os.path.exists(REFERENCE_MAPPING_DB):
        return info
    conn = sqlite3.connect(REFERENCE_MAPPING_DB)
    try:
        c = conn.cursor()
        for gene_id in gene_ids:
            c.execute(
                "SELECT symbol, genename FROM gene_mapping WHERE raw_id = ? LIMIT 1",
                (gene_id,),
            )
            row = c.fetchone()
            if row:
                info[gene_id] = (row[0], row[1])
    finally:
        conn.close()
    return info


def job_genome_for_output(output_dir, default="hg38"):
    """Genome recorded in job.json for the job owning this output directory.

    The parsers need it for two reasons: TargetScan emits one row per species in
    its alignment and only the requested species' rows are targets, and the
    ENST->RefSeq map is per assembly. Guessing human would silently return an
    empty TargetScan result for a mouse job.

    output_dir is <job>/output for a single-pool run and <job>/output/<pool> for
    a network run, so walk up a few levels looking for job.json -- the same
    shape as the targets.txt lookup in result_db._build_db."""
    d = os.path.abspath(output_dir)
    for _ in range(3):
        d = os.path.dirname(d)
        candidate = os.path.join(d, "job.json")
        if os.path.exists(candidate):
            try:
                with open(candidate, "r") as fh:
                    return json.load(fh).get("genome", default) or default
            except (ValueError, IOError):
                return default
    return default


def _build_db(output_dir, db_path):
    """Parse all prediction results and populate a fresh SQLite database.

    Schema:
        gene_tools(gene_id TEXT, tool TEXT, PRIMARY KEY (gene_id, tool))
        gene_info(gene_id TEXT PRIMARY KEY, gene_label TEXT, gene_name TEXT)

    All six tools (including DMISO) are parsed via
    parse_result.process_sequence(). DMISO results are filtered to perfect
    miRNA seed (positions 2-8) matches inside parseDMISOResults().

    Writes to a temp file first then renames atomically so a partial write
    is never visible at db_path.
    """
    json_file = os.path.join(output_dir, "mirna_prediction_parameters.json")
    sequences = read_sequences_from_json(json_file)

    # Discover targets.txt and build ENST->RefSeq map so TargetScan results are
    # converted to RefSeq IDs and, if targets are specified, filtered to the
    # user's list.
    #
    # Prefer one inside output_dir, then fall back to the sibling (job dir).
    # Single-pool workflows write only the sibling, so they are unaffected; a
    # restrict_to_pairs network job writes one per pool inside the pool dir,
    # because each pool has its own target list and there is no single job-level
    # one to share.
    targets_file = os.path.join(output_dir, "targets.txt")
    if not os.path.exists(targets_file):
        targets_file = os.path.join(os.path.dirname(os.path.abspath(output_dir)), "targets.txt")
    targets = load_targets_file(targets_file)
    genome = job_genome_for_output(output_dir)
    enst_to_refseq = build_enst_to_refseq_map(REFERENCE_MAPPING_DB, genome)

    # Collect gene -> tools mapping across all sequences
    gene_tools = {}   # gene_id -> set of tool names

    def _add(gene_id, tool):
        if gene_id is None:
            return
        if gene_id not in gene_tools:
            gene_tools[gene_id] = set()
        gene_tools[gene_id].add(tool)

    for sequence in sequences:
        # Non-DMISO tools via parse_result
        results = process_sequence(
            sequence, output_dir,
            enst_to_refseq=enst_to_refseq, targets=targets,
            genome=genome,
        )
        for tool, gene_ids in results.get("prediction", {}).items():
            for gene_id in gene_ids:
                _add(gene_id, tool)

    # Look up gene labels and names from reference_mapping.db
    gene_info = _load_gene_info(list(gene_tools.keys()))

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
            c.execute("""
                CREATE TABLE gene_info (
                    gene_id    TEXT PRIMARY KEY,
                    gene_label TEXT,
                    gene_name  TEXT
                )
            """)
            c.executemany(
                "INSERT OR IGNORE INTO gene_tools (gene_id, tool) VALUES (?, ?)",
                ((gene_id, tool)
                 for gene_id, tools in gene_tools.items()
                 for tool in tools),
            )
            c.executemany(
                "INSERT OR IGNORE INTO gene_info (gene_id, gene_label, gene_name) VALUES (?, ?, ?)",
                (
                    (gene_id,
                     gene_info.get(gene_id, (None, None))[0],
                     gene_info.get(gene_id, (None, None))[1])
                    for gene_id in gene_tools.keys()
                ),
            )
            conn.commit()
        finally:
            conn.close()
        os.rename(tmp_path, db_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


def _venn_stats(cursor, keyword_pattern=None):
    """Compute per-tool set sizes and all pairwise/higher-order intersections.

    Intersection semantics: |A∩B| counts genes found by BOTH A and B,
    regardless of other tools — the standard input for Venn diagram libraries.

    keyword_pattern: SQL LIKE pattern (e.g. '%TP53%') matched against
                     gene_label, gene_name, or tool name (OR logic). None = no filter.
    """
    cursor.execute("SELECT DISTINCT tool FROM gene_tools ORDER BY tool")
    tools = [row[0] for row in cursor.fetchall()]

    if keyword_pattern:
        gene_filter_sql = (
            "AND gene_id IN ("
            "SELECT gt2.gene_id FROM gene_tools gt2 "
            "LEFT JOIN gene_info gi2 ON gt2.gene_id = gi2.gene_id "
            "WHERE gi2.gene_label LIKE ? OR gi2.gene_name LIKE ? OR gt2.tool LIKE ?)"
        )
        gene_filter_arg = [keyword_pattern, keyword_pattern, keyword_pattern]
    else:
        gene_filter_sql = ""
        gene_filter_arg = []

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

    # Exclusive combinations for the UpSet plot (4+ tools): group genes by the
    # *exact* set of tools that predicted them (so a gene hit by A,B,C lands in
    # one region, not in A&B, A&C, ...). Also roll up per-degree counts -- the
    # number of targets predicted by exactly k tools -- for the consensus view.
    cursor.execute(
        "SELECT gene_id, GROUP_CONCAT(tool) FROM gene_tools "
        "WHERE 1=1 {} GROUP BY gene_id".format(gene_filter_sql),
        gene_filter_arg,
    )
    combo_counts = {}
    for _gene_id, tools_csv in cursor.fetchall():
        if not tools_csv:
            continue
        key = tuple(sorted(set(tools_csv.split(","))))
        combo_counts[key] = combo_counts.get(key, 0) + 1

    combos = []
    degrees = {}
    for key, count in combo_counts.items():
        combos.append({"tools": list(key), "size": count})
        deg = str(len(key))
        degrees[deg] = degrees.get(deg, 0) + count
    combos.sort(key=lambda c: c["size"], reverse=True)

    return {
        "sets": sets,
        "intersections": intersections,
        "combinations": combos,
        "degrees": degrees,
    }


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
                keyword=None):
    """Query the gene-tool database and return paginated results + Venn stats.

    Args:
        sort_by:  "tool_count" | "gene_label"
        order:    "asc" | "desc"
        offset:   row offset for pagination
        number:   page size
        keyword:  optional substring matched (case-insensitive) against gene_label,
                  gene_name, or tool name (OR logic)

    Returns a dict with keys: total_genes, total, genes, venn.
    Each gene entry contains: gene_id, gene_label, gene_name, tool_count, tools.
    """
    sort_dir = "ASC" if order.lower() in ("asc", "ascend") else "DESC"
    sort_col = "tool_count" if sort_by == "tool_count" else "gi.gene_label"

    filter_pattern = "%{}%".format(keyword) if keyword else None

    if filter_pattern:
        where_sql  = (
            "WHERE (gi.gene_label LIKE ? OR gi.gene_name LIKE ? "
            "OR gt.gene_id IN (SELECT gene_id FROM gene_tools WHERE tool LIKE ?))"
        )
        where_args = [filter_pattern, filter_pattern, filter_pattern]
    else:
        where_sql  = ""
        where_args = []

    conn = sqlite3.connect(db_path)
    try:
        c = conn.cursor()

        # Total distinct genes — unfiltered (always the full result set)
        c.execute("SELECT COUNT(DISTINCT gene_id) FROM gene_tools")
        total_genes = c.fetchone()[0]

        # Total distinct genes — respects geneLabel filter
        c.execute(
            """
            SELECT COUNT(DISTINCT gt.gene_id)
            FROM   gene_tools gt
            LEFT JOIN gene_info gi ON gt.gene_id = gi.gene_id
            {}
            """.format(where_sql),
            where_args,
        )
        total = c.fetchone()[0]

        # Paginated gene rows (respects filter)
        c.execute(
            """
            SELECT   gt.gene_id,
                     gi.gene_label,
                     gi.gene_name,
                     COUNT(*)                   AS tool_count,
                     GROUP_CONCAT(gt.tool, ',') AS tools
            FROM     gene_tools gt
            LEFT JOIN gene_info gi ON gt.gene_id = gi.gene_id
            {where}
            GROUP BY gt.gene_id
            ORDER BY {sort_col} {sort_dir}, gt.gene_id ASC
            LIMIT  ? OFFSET ?
            """.format(where=where_sql, sort_col=sort_col, sort_dir=sort_dir),
            where_args + [number, offset],
        )
        genes = [
            {
                "gene_id":    row[0],
                "gene_label": row[1],
                "gene_name":  row[2],
                "tool_count": row[3],
                "tools":      sorted(row[4].split(",")),
            }
            for row in c.fetchall()
        ]

        venn = _venn_stats(c, keyword_pattern=filter_pattern)

    finally:
        conn.close()

    return {"total_genes": total_genes, "total": total, "genes": genes, "venn": venn}


if __name__ == "__main__":
    ensure_db("out")
    res = query_genes("out")
    print(res)
