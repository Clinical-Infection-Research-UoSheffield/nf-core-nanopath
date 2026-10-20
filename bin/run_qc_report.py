#!/usr/bin/env python

"""Build a single run-level QC report from the per-barcode qc_<barcode>.json records.

results_report.py writes one qc_<barcode>.json per barcode that produced a report. This script
rolls them up into one run sheet so the whole run can be quality-checked at a glance from the
pipeline's point of view -- per barcode: reads, clusters formed, clusters that failed to assemble
(and how big), Racon-unpolished clusters, and the fraction of reads lost. Every barcode present in
the run roster but missing a qc json (0 viable reads -> no report) is shown as a failure row, so a
run that analysed nothing can't look blank.

Outputs: <out>.html (human, self-contained) and <out>.tsv (one row per barcode, for auditing/trends).
Pure standard library so it runs in any python3 container.
"""

import argparse
import datetime
import glob
import html
import json
import os

# Per-barcode RAG thresholds (kept in step with results_report.DROPPED_ALARM_PCT).
DROPPED_ALARM_PCT = 30.0

SUMMARY_COLUMNS = [
    ("barcode", "Barcode"),
    ("specimen", "Specimen"),
    ("status", "Status"),
    ("reads_qc", "Reads (QC-passed)"),
    ("reads_clustered", "Reads clustered"),
    ("clusters_formed", "Clusters"),
    ("clusters_classified", "Classified"),
    ("clusters_failed", "Failed"),
    ("clusters_unpolished", "Unpolished"),
    ("pct_reads_dropped", "% reads dropped"),
    ("top_organism", "Top organism"),
    ("top_abundance_pct", "Top %"),
]


def load_roster(path):
    """Full expected barcode list: lines of 'barcode<TAB>status'. Returns [(barcode, status), ...]."""
    roster = []
    if path and os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                parts = line.split("\t")
                roster.append((parts[0].strip(), (parts[1].strip() if len(parts) > 1 else "")))
    return roster


def load_qc(files):
    """barcode -> qc record, from the collected qc_<barcode>.json files."""
    out = {}
    for f in files:
        try:
            with open(f) as fh:
                rec = json.load(fh)
            bc = rec.get("barcode") or os.path.basename(f)[3:-5]
            out[bc] = rec
        except Exception as exc:  # a single bad file must not sink the whole run report
            print("WARN: could not read {0}: {1}".format(f, exc))
    return out


def rag(rec):
    """Red/amber/green for a barcode row. 'missing' records (no report at all) are red."""
    if rec.get("_missing"):
        return "red"
    if int(rec.get("clusters_formed", 0) or 0) == 0:
        return "red"
    if float(rec.get("pct_reads_dropped", 0) or 0) >= DROPPED_ALARM_PCT:
        return "red"
    if (
        int(rec.get("clusters_failed", 0) or 0) > 0
        or int(rec.get("clusters_unpolished", 0) or 0) > 0
        or float(rec.get("pct_reads_dropped", 0) or 0) > 0
    ):
        return "amber"
    return "green"


def _cell(rec, key):
    if rec.get("_missing"):
        return "NO REPORT" if key == "barcode" else ("—" if key != "specimen" else rec.get("specimen", ""))
    v = rec.get(key, "")
    if key == "pct_reads_dropped" and v != "":
        return "{0:.1f}".format(float(v))
    if key == "top_abundance_pct" and v not in ("", None):
        try:
            return "{0:.0f}".format(float(v))
        except (TypeError, ValueError):
            return html.escape(str(v))
    return html.escape(str(v))


CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:#1f2933;margin:24px;background:#fff}
h1{font-size:22px;margin:0 0 4px} .sub{color:#5c6773;margin:0 0 18px;font-size:13px}
.hdr{display:flex;flex-wrap:wrap;gap:8px 22px;margin:0 0 18px;font-size:13px}
.hdr b{color:#1f2933} .tally span{font-weight:700;padding:2px 8px;border-radius:10px;margin-right:6px}
table{border-collapse:collapse;width:100%;font-size:13px;margin:0 0 6px}
th,td{border:1px solid #e3e8ef;padding:5px 8px;text-align:left} th{background:#f5f7fa}
td.c,th.c{text-align:center} tr.red{background:#fdecec} tr.amber{background:#fff7e6}
.dot{display:inline-block;width:10px;height:10px;border-radius:50%;vertical-align:middle}
.dot.red{background:#dc2626}.dot.amber{background:#f59e0b}.dot.green{background:#16a34a}
.pill-red{background:#fdecec;color:#dc2626}.pill-amber{background:#fff7e6;color:#b45309}.pill-green{background:#eafaf0;color:#16a34a}
details{margin:0 0 10px;border:1px solid #e3e8ef;border-radius:6px;padding:4px 10px}
summary{cursor:pointer;font-size:13px;padding:4px 0} summary .dot{margin-right:8px}
details table{margin:8px 0} .banner{background:#dc2626;color:#fff;font-weight:700;padding:10px 14px;border-radius:6px;margin:8px 0}
.foot{color:#8795a5;font-size:11px;margin-top:20px}
"""


def summary_table(rows):
    head = "".join(
        '<th class="{0}">{1}</th>'.format("c" if k not in ("barcode", "specimen", "top_organism") else "", lbl)
        for k, lbl in SUMMARY_COLUMNS
    )
    body = []
    for level, rec in rows:
        tds = ['<td><span class="dot {0}"></span></td>'.format(level)]
        for k, _ in SUMMARY_COLUMNS:
            cls = ' class="c"' if k not in ("barcode", "specimen", "top_organism") else ""
            tds.append("<td{0}>{1}</td>".format(cls, _cell(rec, k)))
        trcls = ' class="{0}"'.format(level) if level in ("red", "amber") else ""
        body.append("<tr{0}>{1}</tr>".format(trcls, "".join(tds)))
    return (
        '<table><thead><tr><th class="c"></th>' + head + "</tr></thead><tbody>" + "".join(body) + "</tbody></table>"
    )


def cluster_details(rec):
    if rec.get("_missing"):
        return ""
    clusters = rec.get("clusters") or []
    if not clusters:
        return ""
    rows = []
    for c in sorted(clusters, key=lambda x: -int(x.get("reads", 0) or 0)):
        failed = c.get("fate") != "classified"
        unpol = bool(c.get("unpolished"))
        trcls = ' class="red"' if failed else (' class="amber"' if unpol else "")
        fate = "unassembled" if failed else ("classified (unpolished)" if unpol else "classified")
        rows.append(
            "<tr{0}><td>{1}</td><td class='c'>{2}</td><td class='c'>{3:.1f}</td>"
            "<td>{4}</td><td>{5}</td></tr>".format(
                trcls,
                html.escape(str(c.get("cluster", ""))),
                int(c.get("reads", 0) or 0),
                float(c.get("pct", 0) or 0),
                fate,
                html.escape(str(c.get("organism", "") or "")),
            )
        )
    level = rag(rec)
    summary = (
        '<span class="dot {0}"></span><b>{1}</b> — {2} cluster(s): '
        "{3} classified, {4} failed, {5} unpolished · {6:.1f}% reads dropped".format(
            level,
            html.escape(str(rec.get("barcode", ""))),
            rec.get("clusters_formed", 0),
            rec.get("clusters_classified", 0),
            rec.get("clusters_failed", 0),
            rec.get("clusters_unpolished", 0),
            float(rec.get("pct_reads_dropped", 0) or 0),
        )
    )
    table = (
        "<table><thead><tr><th>Cluster</th><th class='c'>Reads</th><th class='c'>%</th>"
        "<th>Fate</th><th>Organism</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>"
    )
    return "<details><summary>{0}</summary>{1}</details>".format(summary, table)


def write_tsv(rows, path):
    cols = ["flag"] + [k for k, _ in SUMMARY_COLUMNS]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for level, rec in rows:
            vals = [level]
            for k, _ in SUMMARY_COLUMNS:
                if rec.get("_missing"):
                    vals.append("NO_REPORT" if k == "barcode" else str(rec.get(k, "")))
                else:
                    vals.append(str(rec.get(k, "")))
            fh.write("\t".join(vals) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Build the run-level QC report.")
    ap.add_argument("--qc", nargs="*", default=[], help="per-barcode qc_<barcode>.json files")
    ap.add_argument("--roster", default="", help="TSV of 'barcode<TAB>status' for every expected barcode")
    ap.add_argument("--run_id", default="unknown")
    ap.add_argument("--kit", default="unknown")
    ap.add_argument("--seq_start", default="unknown")
    ap.add_argument("--version", default="unknown")
    ap.add_argument("--out", default="run_qc_report", help="output prefix (.html/.tsv)")
    args = ap.parse_args()

    qc_files = list(args.qc) or sorted(glob.glob("qc_*.json"))
    qc = load_qc(qc_files)
    roster = load_roster(args.roster)
    # If no roster was provided, fall back to whatever qc jsons we have.
    if not roster:
        roster = [(bc, qc[bc].get("status", "")) for bc in sorted(qc)]

    rows = []
    for bc, status in sorted(roster):
        rec = qc.get(bc)
        if rec is None:
            rec = {"barcode": bc, "status": status, "_missing": True}
        rows.append((rag(rec), rec))

    tally = {"red": 0, "amber": 0, "green": 0}
    for level, _ in rows:
        tally[level] = tally.get(level, 0) + 1

    banner = ""
    if tally["red"]:
        banner = (
            '<div class="banner">&#9888; {0} barcode(s) flagged RED — no report produced, no clusters, '
            "or &ge;{1:.0f}% of reads unassembled. Review before sign-off.</div>".format(
                tally["red"], DROPPED_ALARM_PCT
            )
        )

    hdr = (
        '<div class="hdr">'
        "<span><b>Run ID:</b> {run}</span><span><b>Kit:</b> {kit}</span>"
        "<span><b>Seq start:</b> {start}</span><span><b>Barcodes:</b> {n}</span>"
        "<span><b>Pipeline:</b> {ver}</span>"
        "<span><b>Generated:</b> {now}</span>"
        "</div>"
        '<div class="hdr tally"><span class="pill-green">{g} green</span>'
        '<span class="pill-amber">{a} amber</span><span class="pill-red">{r} red</span></div>'
    ).format(
        run=html.escape(args.run_id),
        kit=html.escape(args.kit),
        start=html.escape(args.seq_start),
        n=len(rows),
        ver=html.escape(args.version),
        now=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        g=tally["green"],
        a=tally["amber"],
        r=tally["red"],
    )

    details = "".join(cluster_details(rec) for _, rec in rows)
    doc = (
        "<!doctype html><html><head><meta charset='utf-8'><title>Run QC — {run}</title>"
        "<style>{css}</style></head><body>"
        "<h1>Run QC report</h1>"
        '<p class="sub">Pipeline quality-control summary for the whole run — one row per barcode, '
        "with an expandable per-cluster breakdown. This is <b>not</b> a clinical report.</p>"
        "{hdr}{banner}"
        "<h3>Per-barcode summary</h3>{summary}"
        "<h3>Per-cluster breakdown</h3>{details}"
        '<p class="foot">NanoPATH run QC · failed = cluster formed but produced no consensus (dropped); '
        "unpolished = Racon found too few overlaps and kept the raw draft (lower confidence).</p>"
        "</body></html>"
    ).format(
        run=html.escape(args.run_id),
        css=CSS,
        hdr=hdr,
        banner=banner,
        summary=summary_table(rows),
        details=details or "<p>No per-cluster data.</p>",
    )

    with open(args.out + ".html", "w") as fh:
        fh.write(doc)
    write_tsv(rows, args.out + ".tsv")
    print(
        "run QC report: {0} barcodes ({1} green / {2} amber / {3} red) -> {4}.html, {4}.tsv".format(
            len(rows), tally["green"], tally["amber"], tally["red"], args.out
        )
    )


if __name__ == "__main__":
    main()
