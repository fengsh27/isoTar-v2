"""Did a job actually run the variants it asked for?

mirna_processing.py used to skip a modification whose original base did not
match the sequence (and drop a variant outside 17-30 nt) with only a printed
warning, so a job could run the wild type alone and still report success.
Submission now refuses such a variant and the script runs with --strict, but
jobs finished before that still exist. variant_warnings() compares what a job
requested with the records its mirna.fa really holds and describes anything
missing, so the job page can say the results are wild type only.

mirna.fa headers look like:
    <mirna>,WT
    <mirna>,<mod>[&<mod>...],modified            e.g. hsa-miR-21-5p,3:A|G,modified
    <mirna>,<left>|<right>,shifted               e.g. hsa-let-7a-5p,-1|0,shifted
    <mirna>,<mod>[&<mod>...],<left>|<right>,modified_shifted
"""

VARIANT_KINDS = ("modified", "shifted", "modified_shifted")


def _norm_shift(shift):
    """'+1|0', '1|0' and '1 | 0' are the same shift."""
    try:
        left, right = str(shift).split("|")
        return int(left), int(right)
    except (TypeError, ValueError):
        return str(shift)


def requested_variants(meta):
    """{mirna_id: [{"modifications": [...], "shift": str|None}, ...]} for a job.

    Empty for a custom-sequence job (operations are ignored there) and for a
    job that asked for no variant."""
    if meta.get("workflow") == "mir-network":
        out = {}
        variants = meta.get("mirna_variants")
        if variants is None:
            variants = {}
            for mid in meta.get("mirna_ids") or []:
                specs = []
                s = (meta.get("shifts") or {}).get(mid)
                m = (meta.get("modifications") or {}).get(mid)
                if s or m:
                    specs.append({"shift": s or None, "modifications": list(m or [])})
                specs.extend((meta.get("variants") or {}).get(mid, []))
                if specs:
                    variants[mid] = specs
        for mid, specs in variants.items():
            norm = [{"modifications": list(s.get("modifications") or []),
                     "shift": s.get("shift") or None} for s in specs
                    if s.get("modifications") or s.get("shift")]
            if norm:
                out[mid] = norm
        return out
    if meta.get("mirna_seq"):
        return {}
    mods = list(meta.get("modifications") or [])
    shift = meta.get("shift") or None
    if not (mods or shift):
        return {}
    return {meta.get("mirna_id"): [{"modifications": mods, "shift": shift}]}


def applied_variants(headers):
    """{mirna_id: {"modifications": set, "shifts": set}} from mirna.fa headers."""
    out = {}
    for header in headers:
        header = header.lstrip(">").strip()
        parts = header.split(",")
        if len(parts) < 3 or parts[-1] not in VARIANT_KINDS:
            continue
        entry = out.setdefault(parts[0], {"modifications": set(), "shifts": set()})
        for token in parts[1:-1]:
            if ":" in token:            # modifications ("4:C|U") contain "|" too
                entry["modifications"].update(t for t in token.split("&") if t)
            elif "|" in token:
                entry["shifts"].add(_norm_shift(token))
    return out


def _describe(spec):
    bits = []
    if spec["modifications"]:
        bits.append("modification{} {}".format(
            "s" if len(spec["modifications"]) > 1 else "", ", ".join(spec["modifications"])))
    if spec["shift"]:
        bits.append("shift {}".format(spec["shift"]))
    return " and ".join(bits)


def variant_warnings(meta, headers):
    """Warnings for requested variants missing from the job's mirna.fa.

    Each warning: {code: "variant_not_applied", mirna_id, requested, missing,
    message}. Empty when every requested variant was built."""
    applied = applied_variants(headers)
    warnings = []
    for mid, specs in requested_variants(meta).items():
        got = applied.get(mid, {"modifications": set(), "shifts": set()})
        for spec in specs:
            missing_mods = [m for m in spec["modifications"] if m not in got["modifications"]]
            missing_shift = bool(spec["shift"]) and _norm_shift(spec["shift"]) not in got["shifts"]
            if not missing_mods and not missing_shift:
                continue
            nothing_applied = (len(missing_mods) == len(spec["modifications"])
                               and (missing_shift or not spec["shift"]))
            missing = {"modifications": missing_mods,
                       "shift": spec["shift"] if missing_shift else None}
            if nothing_applied:
                message = ("The requested {} of {} was not applied, so this variant was never "
                           "run: these results are for the wild type only. The original base "
                           "did not match the sequence at that position (or the variant fell "
                           "outside 17-30 nt). Submit the job again to run the variant.").format(
                               _describe(spec), mid)
            else:
                message = ("Part of the requested variant of {} was not applied: {} {} "
                           "missing from the sequence that ran. Submit the job again to run "
                           "the variant as requested.").format(
                               mid, _describe(missing), "is" if len(missing_mods) + missing_shift == 1 else "are")
            warnings.append({"code": "variant_not_applied", "mirna_id": mid,
                             "requested": spec, "missing": missing, "message": message})
    return warnings
