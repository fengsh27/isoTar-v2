"""Identity of the miRNA records a job runs: WT versus engineered variants.

mirna_processing.py writes one FASTA record per sequence it generated, with a
header of the form "<id>,<variant...>":

    hsa-miR-21-5p,WT
    hsa-miR-21-5p,-7|1,shifted
    hsa-miR-21-5p,8:A|U,modified
    hsa-miR-21-5p,8:A|U,-7|1,modified_shifted

Every tool runs once per record and names its raw output after the header, so
the raw results are per sequence. These helpers give each record a stable id
that the result databases (result.db, lncrna_result.db, network.db) store next
to every prediction, so WT and variant hits are never merged.
"""


def mirna_node_id(header):
    """Id for a miRNA record, from the runner header "<id>,<variant...>".

    A WT record keeps the plain base id (so a miRNA with no shift/modification
    is a single node, exactly as before). A modified/shifted variant keeps its
    full header, so it stays distinct from its WT baseline.
        "hsa-miR-21-5p,WT"            -> "hsa-miR-21-5p"
        "hsa-miR-21-5p,-7|1,shifted"  -> "hsa-miR-21-5p,-7|1,shifted"
    """
    parts = header.split(",")
    if len(parts) <= 1 or parts[-1] == "WT":
        return parts[0]
    return header


def mirna_base_id(node_id):
    """Base miRNA id a (possibly variant) node belongs to, for grouping."""
    return node_id.split(",")[0]


def is_variant(node_id):
    """True for a shifted/modified record, False for WT or a plain id."""
    return mirna_base_id(node_id) != node_id


def variant_description(node_id):
    """'WT' for a baseline record, else e.g. 'shifted -7|1' or
    'modified+shifted 8:A|U,-7|1'."""
    parts = node_id.split(",")
    if len(parts) <= 1:
        return "WT"
    vtype = parts[-1].replace("_", "+")
    detail = ",".join(parts[1:-1])
    return "{} {}".format(vtype, detail).strip()


def mirna_node_label(node_id):
    """Human-readable label, e.g. 'hsa-miR-21-5p (shifted -7|1)'. WT/plain nodes
    keep the bare id."""
    if not is_variant(node_id):
        return mirna_base_id(node_id)
    return "{} ({})".format(mirna_base_id(node_id), variant_description(node_id))


def sequence_record(sequence):
    """Row for the `sequences` table from one entry of
    mirna_prediction_parameters.json ({"header", "sequence", ...})."""
    header = sequence["header"]
    node_id = mirna_node_id(header)
    return {
        "sequence": node_id,
        "base_id": mirna_base_id(node_id),
        "header": header,
        "label": mirna_node_label(node_id),
        "variant": variant_description(node_id),
        "is_variant": 1 if is_variant(node_id) else 0,
        "seq": sequence.get("sequence") or "",
    }
