#!/usr/bin/env python
"""Level-0 smoke test for the hit-details feature.

Exercises the pure logic of bin/results_report.py and bin/get_abundance.py against
hand-written fixtures, with aplanat/bokeh/requests stubbed out. No Nextflow, no
containers, no real data are required -- only pandas.

    python -m venv venv && ./venv/bin/pip install pandas
    ./venv/bin/python tests/smoke_hit_details.py        # bin/ is found relative to this file
    ./venv/bin/python tests/smoke_hit_details.py /path/to/bin   # or point at it explicitly

Exits non-zero (via assert) on the first failing check, so it is safe to run in CI.
"""

import sys, os, types, tempfile, importlib.util, csv

# ---- STEP 1: stub the heavy imports and load the two scripts as modules ------
for name in ("aplanat", "aplanat.report", "bokeh", "bokeh.resources", "requests"):
    sys.modules[name] = types.ModuleType(name)
sys.modules["bokeh.resources"].INLINE = None
sys.modules["aplanat"].report = sys.modules["aplanat.report"]

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.join(REPO_ROOT, "bin")


def load(mod_name, filename):
    spec = importlib.util.spec_from_file_location(mod_name, os.path.join(BIN, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


rr = load("results_report", "results_report.py")
ga = load("get_abundance", "get_abundance.py")
print(f"loaded scripts from {BIN}")

# ---- STEP 2: write fixture per-cluster files into a temp dir ------------------
work = tempfile.mkdtemp(prefix="smoke_hits_")
hits = os.path.join(work, "hits")
os.makedirs(hits)


def w(name, text):
    with open(os.path.join(hits, name), "w") as fh:
        fh.write(text)


# cluster 3 = the near-tie: all three classifiers ran (full mode)
#   BLAST columns: sscinames;staxids;evalue;length;pident;bitscore
w(
    "bc_3_blastn_consensus_classification.csv",
    "Streptococcus mitis;28037;0.0;1450;99.10;2600\n"
    "Streptococcus oralis;1303;0.0;1450;99.03;2595\n"
    "Streptococcus pneumoniae;1313;0.0;1448;98.90;2580\n"
    "Streptococcus pseudopneumoniae;257758;0.0;1448;98.20;2550\n",
)  # 4 hits -> top 3 kept
#   SeqMatch columns: sciname;taxid;s_ab_score
w(
    "bc_3_seqmatch_consensus_classification.csv",
    "Streptococcus mitis;28037;0.978\nStreptococcus oralis;1303;0.971\nStreptococcus pneumoniae;1313;0.965\n",
)
#   kraken2 report: pct clade taxon rank taxid name  (tab separated) -> reached species
w(
    "bc_3_kraken2_consensus_classification.csv",
    "  0.00\t0\t0\tU\t0\tunclassified\n100.00\t1\t0\tG\t1301\tStreptococcus\n100.00\t1\t1\tS\t28037\tStreptococcus mitis\n",
)

# cluster 5 = kraken2 only reached genus (the conservative LCA case)
w(
    "bc_5_kraken2_consensus_classification.csv",
    "  0.00\t0\t0\tU\t0\tunclassified\n100.00\t1\t1\tG\t1301\tStreptococcus\n",
)

# cluster 7 = BLAST found nothing (the pipeline's fallback line)
w("bc_7_blastn_consensus_classification.csv", "unclassified;0;0\n")

# cluster 9 = kraken2 left it unclassified
w("bc_9_kraken2_consensus_classification.csv", "100.00\t1\t1\tU\t0\tunclassified\n")

# a chosen-classifier map (what get_abundance would emit in full mode)
chosen = os.path.join(work, "chosen.csv")
with open(chosen, "w") as fh:
    fh.write("cluster,classifier\n3,blast\n5,kraken2\n")

print(f"fixtures in {hits}\n")

# ---- STEP 3: the per-classifier parsers --------------------------------------
blast = rr.parse_blast_hits(os.path.join(hits, "bc_3_blastn_consensus_classification.csv"))
assert list(blast["Species"]) == [
    "Streptococcus mitis *",
    "Streptococcus oralis",
    "Streptococcus pneumoniae",
], "BLAST should keep top-3 in order and asterisk rank 1"
assert list(blast["% Identity"]) == ["99.10", "99.03", "98.90"]
assert list(rr.parse_blast_hits(os.path.join(hits, "bc_7_blastn_consensus_classification.csv"))["Species"]) == [
    "unclassified"
]
kr_genus = rr.parse_kraken2_hit(os.path.join(hits, "bc_5_kraken2_consensus_classification.csv"))
assert list(kr_genus["Species (LCA)"]) == ["Streptococcus *"], "kraken2 collapses to genus, still marked selected"
kr_unc = rr.parse_kraken2_hit(os.path.join(hits, "bc_9_kraken2_consensus_classification.csv"))
assert list(kr_unc["Species (LCA)"]) == ["unclassified"], "unclassified must NOT get the selected asterisk"
print("STEP 3 OK  parsers: top-3, asterisks, genus-collapse, unclassified")


# ---- STEP 4: report grouping + winner-only filtering -------------------------
class FakeSection:
    def __init__(self):
        self.md = []
        self.tables = []

    def markdown(self, s):
        self.md.append(" ".join(s.split()))

    def table(self, df, classes=None):
        self.tables.append(list(df.columns))


class FakeReport:
    def __init__(self):
        self.section = FakeSection()

    def add_section(self):
        return self.section


# without the chosen map -> every classifier shown (single-mode fallback behaviour)
rep = FakeReport()
rr.add_hit_details_section(rep, hits, chosen_classifier="none")
labels_all = [m for m in rep.section.md if m.startswith("_")]
assert "_BLAST_" in labels_all and "_SeqMatch_" in labels_all and "_Kraken2 (LCA)_" in labels_all

# with the chosen map -> cluster 3 shows ONLY blast, cluster 5 ONLY kraken2
rep = FakeReport()
rr.add_hit_details_section(rep, hits, chosen_classifier=chosen)
md = rep.section.md


def labels_after(cluster):
    out, on = [], False
    for m in md:
        if m.startswith("**Cluster"):
            on = m == f"**Cluster {cluster}**"
        elif on and m.startswith("_"):
            out.append(m)
    return out


assert labels_after("3") == ["_BLAST_"], labels_after("3")
assert labels_after("5") == ["_Kraken2 (LCA)_"], labels_after("5")
print("STEP 4 OK  report: all-classifiers without map, winner-only with map")

# ---- STEP 5: the consolidated CSV (top-3 each + top kraken2, ALL clusters) ----
out_csv = os.path.join(work, "hit_details_bc.csv")
rr.write_hit_details_csv(hits, out_csv, chosen_classifier=chosen)
with open(out_csv) as fh:
    rows = list(csv.DictReader(fh))
c3 = [r for r in rows if r["cluster"] == "3"]
assert sorted(r["classifier"] for r in c3) == [
    "BLAST",
    "BLAST",
    "BLAST",
    "Kraken2",
    "SeqMatch",
    "SeqMatch",
    "SeqMatch",
], "cluster 3 must carry all 3 classifiers (3+3+1 rows)"
b1 = [r for r in c3 if r["classifier"] == "BLAST" and r["rank"] == "1"][0]
assert b1["top_hit"] == "True" and b1["pct_identity"] == "99.10" and b1["cluster_winner"] == "BLAST"
assert b1["s_ab_score"] == "" and b1["lca_reads_pct"] == "", "metric columns blank for other classifiers"
print("STEP 5 OK  consolidated CSV: all classifiers, top_hit + cluster_winner + per-metric columns")

# ---- STEP 6: get_abundance winner mapping (+ equivalence to a reference) ------
import pandas as pd

header = (
    "id;reads_in_cluster;used_for_consensus;reads_after_corr;draft_id;"
    "kraken2_sciname;taxid;class_level;name;species;genus;family;order;"
    "seqmatch_sciname;taxid;class_level;name;species;genus;family;order;"
    "blast_sciname;taxid;class_level;name;species;genus;family;order"
)


def block(sci, tid, lvl, nm, sp, ge, fa, orr):
    return ";".join([sci, tid, lvl, nm, sp, ge, fa, orr])


joined = os.path.join(work, "bc.nanoclust_out.txt")
with open(joined, "w") as fh:
    fh.write(header + "\n")
    # cluster 3: kraken2 reached species -> kraken2 wins
    fh.write(
        ";".join(
            [
                "3",
                "120",
                "yes",
                "118",
                "d3",
                block(
                    "S. mitis",
                    "28037",
                    "S",
                    "S. mitis",
                    "S. mitis",
                    "Streptococcus",
                    "Streptococcaceae",
                    "Lactobacillales",
                ),
                block(
                    "S. mitis",
                    "28037",
                    "G",
                    "Streptococcus",
                    "",
                    "Streptococcus",
                    "Streptococcaceae",
                    "Lactobacillales",
                ),
                block(
                    "S. mitis",
                    "28037",
                    "S",
                    "S. mitis",
                    "S. mitis",
                    "Streptococcus",
                    "Streptococcaceae",
                    "Lactobacillales",
                ),
            ]
        )
        + "\n"
    )
    # cluster 5: kraken2 only genus, blast has full species taxonomy -> blast wins
    fh.write(
        ";".join(
            [
                "5",
                "80",
                "yes",
                "78",
                "d5",
                block(
                    "Streptococcus",
                    "1301",
                    "G",
                    "Streptococcus",
                    "",
                    "Streptococcus",
                    "Streptococcaceae",
                    "Lactobacillales",
                ),
                block("E. coli", "562", "G", "Escherichia", "", "Escherichia", "", ""),
                block(
                    "E. coli", "562", "S", "E. coli", "E. coli", "Escherichia", "Enterobacteriaceae", "Enterobacterales"
                ),
            ]
        )
        + "\n"
    )

os.chdir(work)
ga.write_chosen_classifier(joined, "bc")
got = dict((r["cluster"], r["classifier"]) for r in csv.DictReader(open("bc_chosen_classifier.csv")))
assert got == {"3": "kraken2", "5": "blast"}, got

# prove our helper matches the ORIGINAL inline rule from choose_classification
raw = pd.read_csv(joined, index_col=False, sep=";").iloc[:, 1:]


def reference(row):
    if row["class_level"] == "S":
        return "kraken2"
    s = {"kraken2": sum(row.notna()[8:12]), "blast": sum(row.notna()[24:]), "seqmatch": sum(row.notna()[16:20])}
    return max(s, key=s.get)


for _, row in raw.iterrows():
    assert ga.choose_row_classifier(row) == reference(row)
print("STEP 6 OK  winner map correct and equivalent to the original choose_classification rule")

print("\nALL LEVEL-0 CHECKS PASSED")
print("inspect outputs under:", work)
