"""Guard against Python 2.7 SyntaxError on non-ASCII source bytes.

Per CLAUDE.md, v2/*.py scripts must be Python 2.7/3.5 compatible. Python 2
parses source as ASCII by default and raises SyntaxError if a non-ASCII byte
appears in a file with no PEP 263 encoding declaration on line 1 or 2
(e.g. `# -*- coding: utf-8 -*-`). This test catches that class of regression.

Run with:
    python3 -m unittest v2.tests.test_source_encoding
    python2.7 -m unittest v2.tests.test_source_encoding
"""

import io
import os
import re
import unittest


_HERE = os.path.dirname(os.path.abspath(__file__))
_V2_ROOT = os.path.dirname(_HERE)

# PEP 263: encoding declaration must appear on line 1 or 2 and match this regex.
_CODING_RE = re.compile(br"coding[:=][ \t]*([-_.a-zA-Z0-9]+)")


def _iter_v2_py_files():
    for dirpath, _dirnames, filenames in os.walk(_V2_ROOT):
        for name in filenames:
            if name.endswith(".py"):
                yield os.path.join(dirpath, name)


def _has_encoding_declaration(path):
    with io.open(path, "rb") as f:
        head = f.readline() + f.readline()
    return _CODING_RE.search(head) is not None


def _first_non_ascii_line(path):
    """Return (lineno, line_text) of the first line containing a non-ASCII byte,
    or None if the file is pure ASCII."""
    with io.open(path, "rb") as f:
        for i, raw in enumerate(f, start=1):
            try:
                raw.decode("ascii")
            except UnicodeDecodeError:
                # Decode with replacement so we can include a readable snippet.
                return i, raw.decode("utf-8", "replace").rstrip("\r\n")
    return None


class V2SourceEncodingTests(unittest.TestCase):
    """All v2/*.py files must either be pure ASCII or declare an encoding.

    Otherwise `python2.7 /opt/v2/<script>.py` fails with:
        SyntaxError: Non-ASCII character '\\xe2' in file ... but no
        encoding declared; see http://python.org/dev/peps/pep-0263/
    """

    def test_no_non_ascii_without_encoding_declaration(self):
        offenders = []
        for path in _iter_v2_py_files():
            hit = _first_non_ascii_line(path)
            if hit is None:
                continue
            if _has_encoding_declaration(path):
                continue
            lineno, text = hit
            rel = os.path.relpath(path, _V2_ROOT)
            offenders.append("v2/{0}:{1}: {2}".format(rel, lineno, text))

        if offenders:
            self.fail(
                "Found Python file(s) under v2/ with non-ASCII bytes and no "
                "PEP 263 encoding declaration. Either replace the non-ASCII "
                "characters with ASCII equivalents (e.g. em-dash -> --) or "
                "add `# -*- coding: utf-8 -*-` on line 1 or 2:\n  "
                + "\n  ".join(offenders)
            )


if __name__ == "__main__":
    unittest.main()
