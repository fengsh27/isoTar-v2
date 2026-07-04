# -*- coding: utf-8 -*-
"""Canonical lncRNA transcript/gene id normalization for the predictor.

Kept identical to app_v1/lncrna_reference.normalize_lncrna_id so that a target
validated by the backend (POST /api/v1/targets/validate) matches the same record
here when filter_lncrna_fasta subsets the reference FASTA. The two copies exist
because app_v1 (Flask, py3) and v2 (predictor, py2.7/3.x) are separate import
trees -- test_lncrna_validate.py asserts they stay in lockstep. Py2.7-safe.
"""

import re

# Ensembl transcript/gene (ENST/ENSG, ENSMUST/ENSMUSG, ENSCAFT/ENSCAFG, ...) and
# FlyBase (FBtr/FBgn) ids carry a numeric '.<version>' suffix that users usually
# omit. Strip it so 'ENST00000761542' matches the stored 'ENST00000761542.1'.
# WormBase names (Y51H4A.27, WBGene...) have no version and are left intact.
_LNCRNA_ID_RE = re.compile(r'^(ENS[A-Z]{0,6}[TG]\d+|FB(?:tr|gn)\d+)(?:\.\d+)?')


def normalize_lncrna_id(token):
    """Return a lncRNA transcript OR gene id with any Ensembl/FlyBase version
    suffix removed. WormBase names are returned unchanged. None for empty input."""
    if not token:
        return None
    t = token.strip()
    if not t:
        return None
    m = _LNCRNA_ID_RE.match(t)
    if m:
        return m.group(1)
    return t.split()[0]
