#!/usr/bin/env python
"""Regression checks for two linked fixes:

  * abundance percentages must be of the TOTAL reads across ALL clusters the sample formed
    (including dropped/unassembled ones), not just the classified survivors. Otherwise dropping a
    dominant cluster inflates the minor clusters (a ~1% cluster read as 39%).
  * when the unidentified (unassembled) fraction is large (>= DROPPED_ALARM_PCT), the report raises
    a prominent red banner, because a major/dominant organism is then missing from the report.

Needs pandas (./venv/bin/python tests/abundance_denominator.py).
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

# ---- 1) percentages recomputed against the total reads (incl. dropped clusters) ----------------
cluster_info = {
    "3": {"classifier": "BLAST", "reads": 100, "rel_abundance": 50.0},  # inflated: 100/200 survivors
    "5": {"classifier": "BLAST", "reads": 100, "rel_abundance": 50.0},
}
all_clusters = {"1": 9800, "3": 100, "5": 100}  # cluster 1 (dominant) was dropped -> total 10000
rr._recompute_abundance(cluster_info, all_clusters)
assert cluster_info["3"]["rel_abundance"] == 1.0, cluster_info["3"]["rel_abundance"]
assert cluster_info["5"]["rel_abundance"] == 1.0
print("OK  minor cluster 50% (inflated) -> 1.0% of total reads")

ci2 = {"3": {"reads": 100, "rel_abundance": 50.0}}
rr._recompute_abundance(ci2, {})  # no cluster logs available -> leave as-is
assert ci2["3"]["rel_abundance"] == 50.0
print("OK  no cluster logs -> abundance left unchanged (fallback)")


# ---- 2) red banner only when the dropped fraction is large -------------------------------------
def write_logs(mapping):
    d = tempfile.mkdtemp()
    for cid, reads in mapping.items():
        with open(os.path.join(d, "%s.log" % cid), "w") as fh:
            fh.write("%s;%s" % (cid, reads))
    return d


html_big = rr.build_dropped_clusters_html(write_logs({"1": 9800, "3": 100, "5": 100}), ["3", "5"])
assert "background:#dc2626" in html_big and "NOT identified" in html_big and "98%" in html_big, html_big
print("OK  98% unidentified -> red banner raised")

html_small = rr.build_dropped_clusters_html(write_logs({"1": 100, "3": 4900, "5": 5000}), ["3", "5"])
assert "background:#dc2626" not in html_small, "a 1% drop must not raise the red banner"
print("OK  1% unidentified -> note only, no red banner")

print("\nALL ABUNDANCE/BANNER CHECKS PASSED")
