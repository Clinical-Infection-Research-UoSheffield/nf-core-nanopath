#!/usr/bin/env python
"""Checks for the run-level QC report.

Two halves:
  * results_report.write_run_qc_json -- the per-barcode QC record emitted alongside each patient
    report. Verified from SPLIT_CLUSTERS-style cluster logs: cluster counts, the failed
    (unassembled) clusters and their sizes, and the dropped-read fraction. (No classifier fixtures
    needed: with no hit_details, every cluster is 'unassembled', which exercises the counting.)
  * run_qc_report.py -- the aggregator. Driven with synthetic per-barcode qc JSONs plus a roster
    that includes a barcode with NO json, to confirm the RAG logic and that a missing barcode
    surfaces as a red 'NO REPORT' row rather than vanishing.

Needs pandas for the results_report import (./venv/bin/python tests/run_qc_report.py).
"""

import sys, os, json, types, tempfile, importlib.util, contextlib, io

for n in ("aplanat", "aplanat.report", "bokeh", "bokeh.resources", "requests"):
    sys.modules[n] = types.ModuleType(n)
sys.modules["bokeh.resources"].INLINE = None
sys.modules["aplanat"].report = sys.modules["aplanat.report"]
import pandas  # noqa: F401

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "bin")


def _load(mod, fname):
    spec = importlib.util.spec_from_file_location(mod, os.path.join(BIN, fname))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


rr = _load("results_report", "results_report.py")
rq = _load("run_qc_report", "run_qc_report.py")

# ============================ write_run_qc_json ============================
work = tempfile.mkdtemp()
logs = os.path.join(work, "logs")
os.makedirs(logs)
open(os.path.join(logs, "0.log"), "w").write("0;1000")  # big cluster
open(os.path.join(logs, "1.log"), "w").write("1;50")  # small cluster
open(os.path.join(logs, "notacluster.log"), "w").write("ignore me")  # must be skipped

args = types.SimpleNamespace(
    cluster_logs=logs,
    hit_details="none",  # -> no classified clusters, so both clusters are 'unassembled'
    chosen_classifier="none",
    racon_failed="",
    barcode="barcode07",
    status="sample",
    reads_count="1050",
    clustering_size="50000",
    infile="input.1",
)
cwd = os.getcwd()
os.chdir(work)
try:
    rr.write_run_qc_json(args, None)  # metadata_table None is tolerated
    rec = json.load(open("qc_barcode07.json"))
finally:
    os.chdir(cwd)

assert rec["barcode"] == "barcode07", rec
assert rec["clusters_formed"] == 2, rec
assert rec["clusters_classified"] == 0, rec
assert rec["clusters_failed"] == 2, rec
assert rec["reads_clustered"] == 1050, rec
assert rec["pct_reads_dropped"] == 100.0, rec
assert len(rec["clusters"]) == 2 and all(c["fate"] == "unassembled" for c in rec["clusters"]), rec
# biggest cluster first in our emitter's sort key is by id; just confirm sizes captured
sizes = {c["cluster"]: c["reads"] for c in rec["clusters"]}
assert sizes == {"0": 1000, "1": 50}, sizes
print("OK  write_run_qc_json: 2 clusters, both unassembled, sizes + 100% dropped captured")

# ============================ aggregator RAG + missing barcode ============================
qcdir = tempfile.mkdtemp()


def write_qc(rec):
    open(os.path.join(qcdir, "qc_{0}.json".format(rec["barcode"])), "w").write(json.dumps(rec))


def mk(bc, formed, classified, failed, unpol, pct, status="sample"):
    return {
        "barcode": bc,
        "specimen": "S" + bc[-2:],
        "status": status,
        "reads_qc": 5000,
        "reads_clustered": 4800,
        "clustering_cap": 50000,
        "clusters_formed": formed,
        "clusters_classified": classified,
        "clusters_failed": failed,
        "clusters_unpolished": unpol,
        "reads_dropped": int(4800 * pct / 100),
        "pct_reads_dropped": pct,
        "top_organism": "Escherichia coli" if classified else "",
        "top_abundance_pct": 88 if classified else "",
        "clusters": [{"cluster": "0", "reads": 4800, "pct": 100.0, "fate": "classified", "organism": "E. coli",
                      "unpolished": bool(unpol)}] if classified else [],
    }


write_qc(mk("barcode01", formed=2, classified=2, failed=0, unpol=0, pct=0.0))  # green
write_qc(mk("barcode03", formed=5, classified=4, failed=1, unpol=0, pct=5.0))  # amber (a failed cluster)
write_qc(mk("barcode06", formed=40, classified=20, failed=20, unpol=1, pct=42.0))  # red (>30% dropped)

roster = os.path.join(qcdir, "roster.tsv")
with open(roster, "w") as fh:
    fh.write("barcode01\tsample\n")
    fh.write("barcode03\tsample\n")
    fh.write("barcode06\tsample\n")
    fh.write("barcode09\tsample\n")  # NO qc json -> must appear as red NO REPORT row

# direct RAG unit checks
assert rq.rag(mk("x", 2, 2, 0, 0, 0.0)) == "green"
assert rq.rag(mk("x", 5, 4, 1, 0, 5.0)) == "amber"
assert rq.rag(mk("x", 40, 20, 20, 1, 42.0)) == "red"
assert rq.rag({"_missing": True}) == "red"
assert rq.rag(mk("x", 0, 0, 0, 0, 0.0)) == "red"  # zero clusters formed
print("OK  aggregator RAG: green / amber / red / missing / zero-cluster all correct")

out = os.path.join(qcdir, "run_qc_report")
argv = sys.argv
sys.argv = [
    "run_qc_report.py",
    "--qc",
    os.path.join(qcdir, "qc_barcode01.json"),
    os.path.join(qcdir, "qc_barcode03.json"),
    os.path.join(qcdir, "qc_barcode06.json"),
    "--roster", roster,
    "--run_id", "run-XYZ-123",
    "--kit", "SQK-16S114-24",
    "--out", out,
]
try:
    with contextlib.redirect_stdout(io.StringIO()):
        rq.main()
finally:
    sys.argv = argv

doc = open(out + ".html").read()
assert "run-XYZ-123" in doc, "run id missing from header"
assert "NO REPORT" in doc, "missing barcode09 not surfaced as NO REPORT"
assert "Per-cluster breakdown" in doc and "<details>" in doc, "no expandable per-cluster section"
assert "4 red" not in doc  # sanity: not everything red

tsv = [l.split("\t") for l in open(out + ".tsv").read().splitlines()]
header, data = tsv[0], tsv[1:]
assert header[0] == "flag" and "barcode" in header, header
assert len(data) == 4, data  # all four roster barcodes, including the missing one
flag = {r[header.index("barcode")] if r[0] != "red" or "NO_REPORT" not in r else "barcode09": r[0] for r in data}
by_bc = {}
for r in data:
    bc = r[header.index("barcode")]
    by_bc[bc if bc != "NO_REPORT" else "barcode09"] = r[0]
assert by_bc["barcode01"] == "green", by_bc
assert by_bc["barcode03"] == "amber", by_bc
assert by_bc["barcode06"] == "red", by_bc
assert by_bc["barcode09"] == "red", by_bc
print("OK  aggregator: 4 rows (incl. missing barcode as red NO REPORT), flags green/amber/red/red")

print("\nALL RUN-QC-REPORT CHECKS PASSED")
