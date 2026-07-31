#!/usr/bin/env python
"""Prototype check for consensus_cross_contamination.py (contamination-clustering branch).

Two engines are exercised:
  * the built-in edit-distance fallback, on synthetic ~1500 bp sequences (no minimap2 needed here);
  * the minimap2 path, driven by synthetic PAF lines so it needs no minimap2 install -- this is where
    the identity-over-overlap and minimum-overlap rules are verified, including the offset-overlap
    case the built-in fallback misses.

Run with any python3 (edlib optional).
"""

import os, sys, random, tempfile, importlib.util

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("ccc", os.path.join(REPO_ROOT, "bin", "consensus_cross_contamination.py"))
ccc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ccc)

rng = random.Random(42)


def rand_seq(n=1500):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, n_subs):
    s = list(seq)
    for pos in rng.sample(range(len(s)), n_subs):
        s[pos] = rng.choice([b for b in "ACGT" if b != s[pos]])
    return "".join(s)


ecoli = rand_seq()
kleb = rand_seq()
staph = rand_seq()


def write_fasta(records):
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
    for header, seq in records:
        f.write(">{0}\n{1}\n".format(header, seq))
    f.close()
    return f.name


# ============================ built-in fallback (real sequences) ============================
# scenario 1: neg control shares E. coli with a sample -> FLAG
fasta = write_fasta(
    [
        ("barcode01|sample|0|Escherichia coli", ecoli),
        ("barcode03|sample|0|Klebsiella pneumoniae", kleb),
        ("barcode02|negative control|0|Escherichia coli", mutate(ecoli, 5)),
    ]
)
out = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
recs = ccc.load_records(fasta)
neg_hits, sample_hits = ccc.find_matches(recs, min_identity=0.99)  # built-in path (no minimap2 on this host)
assert len(neg_hits) == 1, neg_hits
n, s, pid, ov = neg_hits[0]
assert n["barcode"] == "barcode02" and s["barcode"] == "barcode01", (n["barcode"], s["barcode"])
assert 0.99 <= pid <= 1.0 and ov >= 300, (pid, ov)
assert not sample_hits
print("OK  built-in: neg-control E. coli flagged vs sample at {0:.1f}% over {1} bp".format(pid * 100, ov))

# scenario 2: neg control unrelated to every sample -> CLEAN
fasta2 = write_fasta(
    [
        ("barcode01|sample|0|Escherichia coli", ecoli),
        ("barcode02|negative control|0|unclassified", staph),
    ]
)
assert ccc.find_matches(ccc.load_records(fasta2), 0.99)[0] == []
print("OK  built-in: neg control unrelated to all samples -> no flag")

# scenario 3: two different patients share E. coli -> sample<->sample carryover
fasta3 = write_fasta(
    [
        ("barcode01|sample|0|Escherichia coli", ecoli),
        ("barcode04|sample|0|Escherichia coli", mutate(ecoli, 2)),
        ("barcode05|sample|0|Klebsiella pneumoniae", kleb),
    ]
)
_, sh3 = ccc.find_matches(ccc.load_records(fasta3), 0.99)
assert len(sh3) == 1
print("OK  built-in: two samples sharing E. coli surfaced as sample<->sample")


# ============================ minimap2 path (synthetic PAF) ============================
def rec(barcode, status):
    return {
        "barcode": barcode,
        "status": status,
        "cluster": "0",
        "species": "",
        "seq": "",
        "name": "{0} {1} (cluster 0)".format(barcode, status),
    }


def paf(q, t, matches, block, ln=1500):
    # PAF: qname qlen qstart qend strand tname tlen tstart tend matches block mapq
    return "\t".join(str(x) for x in [q, ln, 0, block, "+", t, ln, 0, block, matches, block, 60])


recs_p = [rec("barcode01", "sample"), rec("barcode02", "negative control"), rec("barcode03", "sample")]

# a) neg(1)~sample(0) at 99.3% over 1500 bp qualifies; neg(1)~sample(2) at 20% is filtered
neg, _ = ccc.find_matches(recs_p, 0.99, 300, paf_lines=[paf(1, 0, 1490, 1500), paf(1, 2, 300, 1500)])
assert len(neg) == 1 and neg[0][0]["barcode"] == "barcode02" and neg[0][1]["barcode"] == "barcode01", neg
assert neg[0][3] == 1500
print("OK  minimap2/PAF: identity-over-overlap flags the neg<->sample match")

# b) minimum-overlap rule: 100% identity but only 200 bp overlap must NOT flag at min_overlap=300
assert ccc.find_matches(recs_p, 0.99, 300, paf_lines=[paf(1, 0, 200, 200)])[0] == []
# ...but the SAME alignment flags once the overlap bar is lowered
assert len(ccc.find_matches(recs_p, 0.99, 150, paf_lines=[paf(1, 0, 200, 200)])[0]) == 1
print("OK  minimap2/PAF: min-overlap rule filters short high-identity hits")

# c) the offset-partial case: two consensus overlap only over 400 bp, but at 99.5% there ->
#    minimap2 reports identity over that overlap, so it is CAUGHT (the built-in would miss it)
neg_off, _ = ccc.find_matches(recs_p, 0.99, 300, paf_lines=[paf(1, 0, 398, 400)])
assert len(neg_off) == 1 and neg_off[0][3] == 400
print("OK  minimap2/PAF: offset/partial overlap (400 bp @ 99.5%) is caught")

# ---- end-to-end: main() writes a TSV (with overlap column) and prints a summary ----
import io, contextlib, argparse

buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    ccc.main(argparse.Namespace(fasta=fasta, min_identity=0.99, min_overlap=300, preset="ava-ont", out=out))
printed = buf.getvalue()
assert "CONTAMINATION SIGNAL" in printed and "barcode02" in printed, printed
header = open(out).read().splitlines()[0]
assert "overlap_bp" in header, header
print("OK  main() writes a TSV with an overlap column and prints a summary")

print("\nALL CROSS-CONTAMINATION CHECKS PASSED")
