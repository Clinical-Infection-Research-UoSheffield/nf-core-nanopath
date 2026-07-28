#!/usr/bin/env python
"""Regression check: check_samplesheet.py rejects duplicate barcodes up front.

SAMPLESHEET_CHECK is the first process in the pipeline. A barcode multiplexes one
physical sample within a run, so it must appear exactly once; a duplicate used to slip
through here and only surface hours later as an opaque crash while building that
barcode's patient report. This asserts the fail-fast, human-readable error instead.

Needs pandas (./venv/bin/python tests/dup_barcode_check.py).
"""

import sys, os, io, tempfile, importlib.util, contextlib
from pathlib import Path
import pandas  # noqa: F401  (check_samplesheet imports pandas at module load)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO_ROOT, "bin")
spec = importlib.util.spec_from_file_location("check_samplesheet", os.path.join(BIN, "check_samplesheet.py"))
cs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cs)


def write_csv(text):
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, newline="")
    f.write(text)
    f.close()
    return Path(f.name)


out = Path(tempfile.NamedTemporaryFile(suffix=".csv", delete=False).name)

# 1) duplicate barcode -> loud SystemExit(1) naming only the duplicate
dup = write_csv("specimen number,barcode,status\n" "S1,barcode05,sample\nS2,barcode05,sample\nS3,barcode06,sample\n")
buf = io.StringIO()
try:
    with contextlib.redirect_stdout(buf):
        cs.check_samplesheet(dup, out, fastq_dir=None, clinical=True)
    raise AssertionError("expected SystemExit on duplicate barcode")
except SystemExit as e:
    assert e.code == 1, e.code
    printed = buf.getvalue()
    assert "ERROR IN THE SAMPLE SHEET - PLEASE REVIEW" in printed, printed
    assert "barcode05" in printed.split("Duplicate")[1], printed
    assert "barcode06" not in printed.split("Duplicate")[1], printed  # only the dup is named
print("OK  duplicate barcode -> loud early SystemExit(1), names only the dup")

# 2) clean sheet -> no exit, writes output
clean = write_csv("specimen number,barcode,status\nS1,barcode05,sample\nS2,barcode06,sample\n")
with contextlib.redirect_stdout(io.StringIO()):
    cs.check_samplesheet(clean, out, fastq_dir=None, clinical=True)  # must NOT raise
assert out.exists()
print("OK  clean sheet passes and writes output")

print("\nALL DUP-BARCODE CHECKS PASSED")
