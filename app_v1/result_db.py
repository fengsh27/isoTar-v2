import csv
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
from app_v1.sequence_ids import sequence_record

DB_FILENAME = "result.db"
PREDICTIONS_CSV = "predictions.csv"

# Bumped whenever the tables change shape. ensure_db() rebuilds any cached
# database carrying an older `PRAGMA user_version` from the raw tool output,
# which is still on disk, so a schema change needs no migration script.
#   1: gene_tools(gene_id, tool) -- WT and variant hits merged (implicit, no
#      user_version written)
#   2: sequences + gene_tools(sequence, gene_id, tool) -- one row per record
#      that was actually run, so WT and variants can be told apart
SCHEMA_VERSION = 2

# Path to the reference mapping database (gene_id -> gene_label/gene_name)
REFERENCE_MAPPING_DB = os.environ.get(
    "ISOTAR_REFERENCE_MAPPING_DB",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference_mapping.db"),
)

COMPARE_STATUSES = ("changed", "lost", "gained", "shared")


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


def collect_gene_edges(output_dir, enst_to_refseq, targets=None, genome="hg38"):
    """Parse every record's raw tool output in output_dir.

    Returns (sequence_rows, edges):
        sequence_rows: one sequence_record() per record in
                       mirna_prediction_parameters.json, in file order
        edges:         sequence id -> {gene_id -> set(tool names)}

    All six tools (including DMISO) are parsed via
    parse_result.process_sequence(), per record, so a variant's hits stay
    attributed to the variant. DMISO results are filtered to perfect miRNA
    seed (positions 2-8) matches inside parseDMISOResults().
    """
    json_file = os.path.join(output_dir, "mirna_prediction_parameters.json")
    sequence_rows = []
    edges = {}
    for sequence in read_sequences_from_json(json_file):
        row = sequence_record(sequence)
        sequence_rows.append(row)
        targets_of = edges.setdefault(row["sequence"], {})
        results = process_sequence(
            sequence, output_dir,
            enst_to_refseq=enst_to_refseq, targets=targets,
            genome=genome,
        )
        for tool, gene_ids in results.get("prediction", {}).items():
            for gene_id in gene_ids:
                if gene_id:
                    targets_of.setdefault(gene_id, set()).add(tool)
    return sequence_rows, edges


def write_db(output_dir, db_path, sequence_rows, edges, gene_info, csv_name=None):
    """Write a fresh per-sequence result database (and optionally its CSV).

    Schema (SCHEMA_VERSION 2):
        sequences(sequence PK, base_id, header, label, variant, is_variant, seq)
        gene_tools(sequence, gene_id, tool, PRIMARY KEY (sequence, gene_id, tool))
        gene_info(gene_id PK, gene_label, gene_name)

    gene_info: gene_id -> (gene_label, gene_name); targets missing from it get
    NULL label and name.

    Writes to a temp file first then renames atomically so a partial write is
    never visible at db_path, and a reader holding the old file keeps it.
    """
    gene_ids = set()
    for targets_of in edges.values():
        gene_ids.update(targets_of.keys())

    tmp_fd, tmp_path = tempfile.mkstemp(dir=output_dir, suffix=".db.tmp")
    os.close(tmp_fd)
    try:
        conn = sqlite3.connect(tmp_path)
        try:
            c = conn.cursor()
            c.execute("""
                CREATE TABLE sequences (
                    sequence   TEXT PRIMARY KEY,
                    base_id    TEXT NOT NULL,
                    header     TEXT NOT NULL,
                    label      TEXT NOT NULL,
                    variant    TEXT NOT NULL,
                    is_variant INTEGER NOT NULL,
                    seq        TEXT NOT NULL
                )
            """)
            c.execute("""
                CREATE TABLE gene_tools (
                    sequence TEXT NOT NULL,
                    gene_id  TEXT NOT NULL,
                    tool     TEXT NOT NULL,
                    PRIMARY KEY (sequence, gene_id, tool)
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
                "INSERT OR IGNORE INTO sequences "
                "(sequence, base_id, header, label, variant, is_variant, seq) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                ((r["sequence"], r["base_id"], r["header"], r["label"],
                  r["variant"], r["is_variant"], r["seq"]) for r in sequence_rows),
            )
            c.executemany(
                "INSERT OR IGNORE INTO gene_tools (sequence, gene_id, tool) VALUES (?, ?, ?)",
                ((sequence, gene_id, tool)
                 for sequence, targets_of in edges.items()
                 for gene_id, tools in targets_of.items()
                 for tool in tools),
            )
            c.executemany(
                "INSERT OR IGNORE INTO gene_info (gene_id, gene_label, gene_name) VALUES (?, ?, ?)",
                ((gene_id,
                  gene_info.get(gene_id, (None, None))[0],
                  gene_info.get(gene_id, (None, None))[1])
                 for gene_id in gene_ids),
            )
            c.execute("PRAGMA user_version = {}".format(SCHEMA_VERSION))
            conn.commit()
        finally:
            conn.close()
        os.rename(tmp_path, db_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise

    if csv_name:
        _write_csv(output_dir, os.path.join(output_dir, csv_name),
                   sequence_rows, edges, gene_info)


def _write_csv(output_dir, csv_path, sequence_rows, edges, gene_info):
    """Flat export of the same rows, one prediction per line, for users who
    want the browser table outside the browser. Sorted so repeated builds of
    the same job give byte-identical files."""
    by_id = dict((r["sequence"], r) for r in sequence_rows)
    tmp_fd, tmp_path = tempfile.mkstemp(dir=output_dir, suffix=".csv.tmp")
    os.close(tmp_fd)
    try:
        with open(tmp_path, "w") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["sequence", "mirna", "variant", "gene_id",
                        "gene_label", "gene_name", "tool"])
            for sequence in sorted(edges):
                rec = by_id.get(sequence, {})
                for gene_id in sorted(edges[sequence]):
                    label, name = gene_info.get(gene_id, (None, None))
                    for tool in sorted(edges[sequence][gene_id]):
                        w.writerow([sequence, rec.get("base_id", sequence),
                                    rec.get("variant", ""), gene_id,
                                    label or "", name or "", tool])
        os.rename(tmp_path, csv_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


def _build_db(output_dir, db_path):
    """Parse all prediction results and populate a fresh SQLite database.

    See write_db() for the schema. Also writes PREDICTIONS_CSV next to it.
    """
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

    sequence_rows, edges = collect_gene_edges(
        output_dir, enst_to_refseq, targets=targets, genome=genome)

    gene_ids = set()
    for targets_of in edges.values():
        gene_ids.update(targets_of.keys())
    # Look up gene labels and names from reference_mapping.db
    gene_info = _load_gene_info(list(gene_ids))

    write_db(output_dir, db_path, sequence_rows, edges, gene_info,
             csv_name=PREDICTIONS_CSV)


def schema_version(db_path):
    """`PRAGMA user_version` of a cached database; 0 for a version-1 file
    (which never wrote one) and -1 if the file cannot be opened."""
    try:
        conn = sqlite3.connect(db_path)
    except sqlite3.Error:
        return -1
    try:
        return conn.execute("PRAGMA user_version").fetchone()[0]
    except sqlite3.Error:
        return -1
    finally:
        conn.close()


class _SequenceFilter(object):
    """SQL fragment restricting gene_tools rows to a set of sequence ids.

    `sql(alias)` renders "AND <alias>sequence IN (?, ?)" for the given table
    alias prefix ("gt." or ""), or "" when no restriction applies -- a
    version-1 database has no sequence column, so it is always unrestricted.
    """

    def __init__(self, sequences):
        self.sequences = list(sequences) if sequences is not None else None

    @property
    def active(self):
        return self.sequences is not None

    def sql(self, alias=""):
        if not self.active:
            return ""
        if not self.sequences:
            return "AND 0"
        return "AND {}sequence IN ({})".format(alias, ",".join("?" * len(self.sequences)))

    @property
    def args(self):
        return list(self.sequences) if self.sequences else []


def _has_sequences(cursor):
    cursor.execute("SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'sequences'")
    return cursor.fetchone() is not None


def _wt_sequences(cursor):
    cursor.execute("SELECT sequence FROM sequences WHERE is_variant = 0 ORDER BY sequence")
    return [row[0] for row in cursor.fetchall()]


def _resolve_sequences(cursor, sequences):
    """Turn the caller's `sequences` argument into a _SequenceFilter.

    None means "every WT record" -- the baseline a user expects by default,
    and exactly what a job without variants always showed. A job whose only
    records are variants (no WT survived) falls back to all records rather
    than an empty page. A version-1 database is never filtered.
    """
    if not _has_sequences(cursor):
        return _SequenceFilter(None)
    if sequences is None:
        wt = _wt_sequences(cursor)
        if not wt:
            return _SequenceFilter(None)
        return _SequenceFilter(wt)
    return _SequenceFilter(sequences)


def _venn_stats(cursor, keyword_pattern=None, seq_filter=None):
    """Compute per-tool set sizes and all pairwise/higher-order intersections.

    Intersection semantics: |A∩B| counts genes found by BOTH A and B,
    regardless of other tools — the standard input for Venn diagram libraries.

    keyword_pattern: SQL LIKE pattern (e.g. '%TP53%') matched against
                     gene_label, gene_name, or tool name (OR logic). None = no filter.
    seq_filter:      _SequenceFilter restricting which records' hits count.
    """
    seq_filter = seq_filter or _SequenceFilter(None)
    seq_sql = seq_filter.sql()
    seq_args = seq_filter.args

    cursor.execute(
        "SELECT DISTINCT tool FROM gene_tools WHERE 1=1 {} ORDER BY tool".format(seq_sql),
        seq_args,
    )
    tools = [row[0] for row in cursor.fetchall()]

    if keyword_pattern:
        gene_filter_sql = (
            "AND gene_id IN ("
            "SELECT gt2.gene_id FROM gene_tools gt2 "
            "LEFT JOIN gene_info gi2 ON gt2.gene_id = gi2.gene_id "
            "WHERE (gi2.gene_label LIKE ? OR gi2.gene_name LIKE ? OR gt2.tool LIKE ?) {})"
        ).format(seq_filter.sql("gt2."))
        gene_filter_arg = [keyword_pattern, keyword_pattern, keyword_pattern] + seq_args
    else:
        gene_filter_sql = ""
        gene_filter_arg = []

    sets = {}
    for tool in tools:
        cursor.execute(
            "SELECT COUNT(DISTINCT gene_id) FROM gene_tools WHERE tool = ? {} {}".format(
                gene_filter_sql, seq_sql
            ),
            [tool] + gene_filter_arg + seq_args,
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
                    WHERE  tool IN ({}) {} {}
                    GROUP  BY gene_id
                    HAVING COUNT(DISTINCT tool) = ?
                )
                """.format(placeholders, gene_filter_sql, seq_sql),
                list(combo) + gene_filter_arg + seq_args + [len(combo)],
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
        "WHERE 1=1 {} {} GROUP BY gene_id".format(gene_filter_sql, seq_sql),
        gene_filter_arg + seq_args,
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
    """Return path to result.db, building it first if it does not exist or was
    written by an older schema (see SCHEMA_VERSION)."""
    db_path = os.path.join(output_dir, DB_FILENAME)
    if not os.path.exists(db_path) or schema_version(db_path) < SCHEMA_VERSION:
        _build_db(output_dir, db_path)
    return db_path


def list_sequences(db_path):
    """The records the job ran, WT first within each miRNA, for the UI's
    sequence selector. Empty for a version-1 database."""
    conn = sqlite3.connect(db_path)
    try:
        c = conn.cursor()
        if not _has_sequences(c):
            return []
        c.execute(
            "SELECT sequence, base_id, header, label, variant, is_variant, seq "
            "FROM sequences ORDER BY base_id, is_variant, sequence"
        )
        return [
            {
                "sequence":   row[0],
                "base_id":    row[1],
                "header":     row[2],
                "label":      row[3],
                "variant":    row[4],
                "is_variant": bool(row[5]),
                "seq":        row[6],
            }
            for row in c.fetchall()
        ]
    finally:
        conn.close()


def gene_labels(db_path, sequences=None):
    """Every distinct gene symbol among the predicted targets, sorted.

    The input to enrichment: targets without a symbol (gene_info has no label
    for them) are skipped, and transcripts of one gene collapse to its symbol.
    `sequences` picks the records whose targets count; None = the WT records,
    so a variant's hypothetical targets never leak into the enrichment of the
    canonical miRNA.
    """
    conn = sqlite3.connect(db_path)
    try:
        c = conn.cursor()
        seq_filter = _resolve_sequences(c, sequences)
        c.execute(
            """
            SELECT DISTINCT TRIM(gi.gene_label)
            FROM   gene_tools gt
            JOIN   gene_info gi ON gt.gene_id = gi.gene_id
            WHERE  gi.gene_label IS NOT NULL AND TRIM(gi.gene_label) != '' {}
            ORDER BY 1
            """.format(seq_filter.sql("gt.")),
            seq_filter.args,
        )
        return [row[0] for row in c.fetchall()]
    finally:
        conn.close()


def query_genes(db_path, sort_by="tool_count", order="desc", offset=0, number=20,
                keyword=None, sequences=None):
    """Query the gene-tool database and return paginated results + Venn stats.

    Args:
        sort_by:   "tool_count" | "gene_label"
        order:     "asc" | "desc"
        offset:    row offset for pagination
        number:    page size
        keyword:   optional substring matched (case-insensitive) against gene_label,
                   gene_name, or tool name (OR logic)
        sequences: ids of the records whose hits to show (see list_sequences).
                   None = the WT records. Several ids give their union, with a
                   tool counted once per gene however many records hit it.

    Returns a dict with keys: total_genes, total, genes, venn, selected.
    Each gene entry contains: gene_id, gene_label, gene_name, tool_count, tools.
    """
    sort_dir = "ASC" if order.lower() in ("asc", "ascend") else "DESC"
    sort_col = "tool_count" if sort_by == "tool_count" else "gi.gene_label"

    filter_pattern = "%{}%".format(keyword) if keyword else None

    conn = sqlite3.connect(db_path)
    try:
        c = conn.cursor()
        seq_filter = _resolve_sequences(c, sequences)
        seq_sql = seq_filter.sql("gt.")
        seq_args = seq_filter.args

        if filter_pattern:
            where_sql = (
                "WHERE (gi.gene_label LIKE ? OR gi.gene_name LIKE ? "
                "OR gt.gene_id IN (SELECT gt3.gene_id FROM gene_tools gt3 "
                "WHERE gt3.tool LIKE ? {})) {}"
            ).format(seq_filter.sql("gt3."), seq_sql)
            where_args = [filter_pattern, filter_pattern, filter_pattern] + seq_args + seq_args
        else:
            where_sql = "WHERE 1=1 {}".format(seq_sql)
            where_args = list(seq_args)

        # Total distinct genes — unfiltered by keyword (always the full result set)
        c.execute(
            "SELECT COUNT(DISTINCT gt.gene_id) FROM gene_tools gt WHERE 1=1 {}".format(seq_sql),
            seq_args,
        )
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
                     COUNT(DISTINCT gt.tool)       AS tool_count,
                     GROUP_CONCAT(DISTINCT gt.tool) AS tools
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

        venn = _venn_stats(c, keyword_pattern=filter_pattern, seq_filter=seq_filter)
        selected = seq_filter.sequences if seq_filter.active else []

    finally:
        conn.close()

    return {"total_genes": total_genes, "total": total, "genes": genes,
            "venn": venn, "selected": selected}


def compare_genes(db_path, wt, variant, status=None, keyword=None,
                  sort_by="status", order="desc", offset=0, number=20):
    """Per-target comparison of one variant record against its WT record.

    Every gene hit by either record is classified:
        shared   same tools in both
        gained   no WT tools, at least one variant tool  (variant-only target)
        lost     WT tools, no variant tools              (target the variant loses)
        changed  both non-empty but different; `gained`/`lost` on the row list
                 which tools moved
    "Lost" means no longer predicted by those tools, not evidence of no
    binding -- the UI should say so.

    Args:
        status:  restrict rows to one COMPARE_STATUSES value
        keyword: case-insensitive substring of gene_id, gene_label or gene_name
        sort_by: "status" (changed, lost, gained, shared; then most tools moved)
                 | "gene_label" | "delta" (number of tools moved)
        order:   "asc" | "desc" -- "desc" is the natural order for "status"

    Returns: wt, variant, total_genes (all genes of either record), total
    (after filters), totals (genes per status, unfiltered), tools (per tool:
    genes where it was gained / lost), genes (the page).
    """
    conn = sqlite3.connect(db_path)
    try:
        c = conn.cursor()
        c.execute(
            """
            SELECT gt.gene_id, gi.gene_label, gi.gene_name, gt.sequence, gt.tool
            FROM   gene_tools gt
            LEFT JOIN gene_info gi ON gt.gene_id = gi.gene_id
            WHERE  gt.sequence IN (?, ?)
            """,
            (wt, variant),
        )
        rows = {}
        for gene_id, label, name, sequence, tool in c.fetchall():
            r = rows.get(gene_id)
            if r is None:
                r = rows[gene_id] = {
                    "gene_id": gene_id, "gene_label": label, "gene_name": name,
                    "_wt": set(), "_variant": set(),
                }
            (r["_wt"] if sequence == wt else r["_variant"]).add(tool)
    finally:
        conn.close()

    totals = dict((s, 0) for s in COMPARE_STATUSES)
    per_tool = {}
    for r in rows.values():
        wt_tools, var_tools = r.pop("_wt"), r.pop("_variant")
        gained, lost = var_tools - wt_tools, wt_tools - var_tools
        if not wt_tools:
            r["status"] = "gained"
        elif not var_tools:
            r["status"] = "lost"
        elif gained or lost:
            r["status"] = "changed"
        else:
            r["status"] = "shared"
        r["wt_tools"] = sorted(wt_tools)
        r["variant_tools"] = sorted(var_tools)
        r["gained"] = sorted(gained)
        r["lost"] = sorted(lost)
        r["delta"] = len(gained) + len(lost)
        totals[r["status"]] += 1
        for tool in gained:
            per_tool.setdefault(tool, {"gained": 0, "lost": 0})["gained"] += 1
        for tool in lost:
            per_tool.setdefault(tool, {"gained": 0, "lost": 0})["lost"] += 1

    genes = list(rows.values())
    if status:
        genes = [g for g in genes if g["status"] == status]
    if keyword:
        needle = keyword.lower()
        genes = [g for g in genes
                 if needle in (g["gene_id"] or "").lower()
                 or needle in (g["gene_label"] or "").lower()
                 or needle in (g["gene_name"] or "").lower()]

    reverse = order.lower() not in ("asc", "ascend")
    rank = dict((s, len(COMPARE_STATUSES) - i) for i, s in enumerate(COMPARE_STATUSES))
    # Label order is the tiebreak in every mode and always reads A..Z, so sort
    # it first (sorts are stable) and let `order` flip only the primary key.
    genes.sort(key=lambda g: ((g["gene_label"] or g["gene_id"]).lower(), g["gene_id"]))
    if sort_by == "gene_label":
        if reverse:
            genes.reverse()
    elif sort_by == "delta":
        genes.sort(key=lambda g: g["delta"], reverse=reverse)
    else:
        genes.sort(key=lambda g: (rank[g["status"]], g["delta"]), reverse=reverse)

    total = len(genes)
    page = genes[offset:offset + number]

    return {
        "wt": wt,
        "variant": variant,
        "total_genes": len(rows),
        "total": total,
        "totals": totals,
        "tools": per_tool,
        "genes": page,
    }


if __name__ == "__main__":
    ensure_db("out")
    res = query_genes(os.path.join("out", DB_FILENAME))
    print(res)
