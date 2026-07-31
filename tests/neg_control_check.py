#!/usr/bin/env python
"""Regression check: the whole-run negative-control read count.

The negative control is no longer a per-cluster traffic light. It is assessed once for the whole
run: any reads at all in the negative control -> flagged (red) with a note to read the sample
report alongside the negative-control report. negative_control_reads sums reads across the full
neg-control file (not the top-3 display table).

Needs pandas (./venv/bin/python tests/neg_control_check.py).
"""

import sys, os, types, importlib.util, tempfile

for n in ("aplanat", "aplanat.report", "bokeh", "bokeh.resources", "requests"):
    sys.modules[n] = types.ModuleType(n)
sys.modules["bokeh.resources"].INLINE = None
sys.modules["aplanat"].report = sys.modules["aplanat.report"]
import pandas  # noqa: F401

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("rr", os.path.join(REPO_ROOT, "bin", "results_report.py"))
rr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rr)


def control_file(text):
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False)
    f.write(text)
    f.close()
    return f.name


# neg control WITH reads across several species (incl. one beyond the top 3) -> full total counted
withreads = control_file(
    "taxid,rel_abundance,reads\n"
    "Escherichia coli,40,400\nCutibacterium acnes,30,300\nBacillus subtilis,20,200\n"
    "Staphylococcus epidermidis,10,100\n"
)
assert rr.negative_control_reads("[%s]" % withreads) == 1000, rr.negative_control_reads("[%s]" % withreads)
print("OK  sums reads across the full file (400+300+200+100 = 1000), not just the top 3")

# neg control with zero reads -> 0 (green)
zero = control_file("taxid,rel_abundance,reads\n")
assert rr.negative_control_reads("[%s]" % zero) == 0
print("OK  empty negative control -> 0 reads")

# no negative control / failed / [None] -> 0, no crash
assert rr.negative_control_reads("[None]") == 0
assert rr.negative_control_reads("[/no/such/file.csv]") == 0
print("OK  missing / [None] control -> 0 reads, no crash")

print("\nALL NEG-CONTROL CHECKS PASSED")
