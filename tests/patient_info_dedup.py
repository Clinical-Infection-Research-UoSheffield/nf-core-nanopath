#!/usr/bin/env python
"""Regression checks for read_patient_info() duplicate-barcode handling.

A barcode identifies one physical sample, so the samplesheet must have one row per
barcode. The old code transposed *all* matching rows and hard-coded two column names,
so a duplicated barcode crashed with an opaque pandas "Length mismatch: Expected axis
has 3 elements, new values have 2 elements". This asserts the new behaviour:

  * one row              -> normal 2-column Metadata/Sample-Information frame
  * exact-duplicate rows -> collapsed (tolerated: harmless copy/paste)
  * conflicting rows     -> clear, actionable SystemExit (unsafe to guess the patient)
  * missing barcode      -> unchanged not-found error

Needs pandas (./venv/bin/python tests/patient_info_dedup.py). read_excel is stubbed so
no Excel engine (openpyxl) is required.
"""

import sys, os, types, importlib.util

for name in ("aplanat", "aplanat.report", "bokeh", "bokeh.resources", "requests"):
    sys.modules[name] = types.ModuleType(name)
sys.modules["bokeh.resources"].INLINE = None
sys.modules["aplanat"].report = sys.modules["aplanat.report"]
import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO_ROOT, "bin")
spec = importlib.util.spec_from_file_location("results_report", os.path.join(BIN, "results_report.py"))
rr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rr)

COLS = ["Specimen Number", "Barcode", "Status", "Name"] + ["c%d" % i for i in range(7)]


def as_df(rows):
    return pd.DataFrame(rows, columns=COLS)


# feed a DataFrame straight in as "file"; stub read_excel to mimic the real usecols=range(0,11)
rr.pd.read_excel = lambda f, usecols=None, **kw: (f.iloc[:, usecols] if usecols is not None else f)

base = ["S1", "barcode05", "sample", "Patient A"] + [""] * 7

df1 = as_df([base, ["S2", "barcode06", "sample", "Patient B"] + [""] * 7])
out = rr.read_patient_info(df1, "barcode05")
assert list(out.columns) == ["Metadata", "Sample Information"], out.columns
assert "barcode05" in out["Sample Information"].astype(str).tolist()
print("OK  single row -> 2-column frame")

out2 = rr.read_patient_info(as_df([base, list(base)]), "barcode05")
assert list(out2.columns) == ["Metadata", "Sample Information"]
print("OK  exact-duplicate rows collapse cleanly")

try:
    rr.read_patient_info(as_df([base, ["S9", "barcode05", "sample", "Patient DIFFERENT"] + [""] * 7]), "barcode05")
    raise AssertionError("expected SystemExit on conflicting duplicate")
except SystemExit as e:
    msg = str(e)
    assert "appears 2 times" in msg and "barcode05" in msg, msg
    assert "S1" in msg and "S9" in msg, msg
print("OK  conflicting duplicate -> clear, actionable error")

try:
    rr.read_patient_info(df1, "barcode99")
    raise AssertionError("expected SystemExit on missing barcode")
except SystemExit as e:
    assert "not found" in str(e)
print("OK  missing barcode -> not-found error")

print("\nALL PATIENT-INFO CHECKS PASSED")
