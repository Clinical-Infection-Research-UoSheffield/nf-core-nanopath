#!/usr/bin/env python
"""Level-0 checks for the two QC footers in build_qc_html:

  * "unshown clusters" note   -- clusters below the abundance threshold are counted, not dropped
  * "Unpolished consensus"    -- clusters whose Racon fell back to the raw draft are flagged amber

No pandas / Nextflow needed: build_qc_html works on plain dicts, so we stub the heavy
imports (same trick as smoke_hit_details.py) and call it directly.

    python tests/qc_footers.py

Exits non-zero (via assert) on the first failing check.
"""

import sys, os, types, importlib.util

for name in ("aplanat", "aplanat.report", "bokeh", "bokeh.resources", "requests", "pandas", "numpy"):
    sys.modules[name] = types.ModuleType(name)
sys.modules["bokeh.resources"].INLINE = None
sys.modules["aplanat"].report = sys.modules["aplanat.report"]

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO_ROOT, "bin")
spec = importlib.util.spec_from_file_location("results_report", os.path.join(BIN, "results_report.py"))
rr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rr)


def krak(species):
    """A minimal kraken2 rank-1 record; other classifiers absent (that's fine for these checks)."""
    return {"species": species, "pct_identity": "", "bitscore": "", "s_ab_score": ""}


# cluster 0 = big (shown);  cluster 1 = tiny (below the 5% threshold -> unshown)
clusters = {"0": {"kraken2": [krak("Streptococcus mitis")]}, "1": {"kraken2": [krak("Escherichia coli")]}}
cluster_info = {
    "0": {"reads": 900, "rel_abundance": 90.0, "classifier": "kraken2"},
    "1": {"reads": 30, "rel_abundance": 3.0, "classifier": "kraken2"},
}

# ---- 1) unshown footer appears and names the count + largest hidden abundance -----------
html = rr.build_qc_html(clusters, cluster_info, neg_species=[])
assert "additional cluster(s) which are not shown" in html, "unshown-clusters note missing"
assert "The largest was 3.0%" in html, "unshown note should report the largest hidden abundance"
assert "1</td>" in html or ">1<" in html or "There were 1 " in html, "should count the 1 hidden cluster"
# the hidden cluster's species must NOT be rendered as a row
assert "Escherichia coli" not in html, "below-threshold cluster must not be shown in the table"
print("OK  unshown-clusters footer: counts hidden clusters + largest abundance, hides the row")

# ---- 2) no unshown footer when everything clears the threshold --------------------------
html_all = rr.build_qc_html({"0": clusters["0"]}, {"0": cluster_info["0"]}, neg_species=[])
assert "additional cluster(s) which are not shown" not in html_all, "no hidden clusters -> no note"
print("OK  unshown-clusters footer absent when nothing is hidden")

# ---- 3) Racon footer appears only for the flagged cluster -------------------------------
html_racon = rr.build_qc_html(clusters, cluster_info, neg_species=[], racon_failed=frozenset({"0"}))
assert "Unpolished consensus" in html_racon, "racon-failed footer missing"
assert "Cluster(s) 0" in html_racon, "racon footer should name the affected cluster"
print("OK  Racon footer present and names the affected cluster")

# ---- 4) no Racon footer when the set is empty (default) ---------------------------------
assert "Unpolished consensus" not in html, "no racon-failed clusters -> no unpolished-consensus note"
print("OK  Racon footer absent by default")

print("\nALL QC-FOOTER CHECKS PASSED")
