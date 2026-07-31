#!/usr/bin/env python
"""Prototype check for consensus_cross_contamination.py (contamination-clustering branch).

Builds a small synthetic run: a sample and a negative control that share a near-identical ~1500 bp
consensus (planted contamination), plus unrelated sequences that must not flag. Verifies the
sequence-level check flags the neg<->sample match and leaves the clean sequences alone.

Needs edlib (./venv/bin/pip install edlib) and is run with the venv python.
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


ecoli = rand_seq()  # the organism shared between a sample and the neg control
kleb = rand_seq()  # an unrelated organism, genuinely only in a sample
staph = rand_seq()  # unrelated, only in the neg control in the clean scenario


def write_fasta(records):
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
    for header, seq in records:
        f.write(">{0}\n{1}\n".format(header, seq))
    f.close()
    return f.name


# ---- scenario 1: neg control shares E. coli with a sample (5 subs = 99.7%) -> FLAG -------------
fasta = write_fasta(
    [
        ("barcode01|sample|0|Escherichia coli", ecoli),
        ("barcode03|sample|0|Klebsiella pneumoniae", kleb),
        ("barcode02|negative control|0|Escherichia coli", mutate(ecoli, 5)),
    ]
)
out = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
recs = ccc.load_records(fasta)
neg_hits, sample_hits = ccc.find_matches(recs, min_identity=0.99)
assert len(neg_hits) == 1, neg_hits
n, s, pid = neg_hits[0]
assert n["barcode"] == "barcode02" and s["barcode"] == "barcode01", (n["barcode"], s["barcode"])
assert 0.99 <= pid <= 1.0, pid
assert not sample_hits, "unrelated sample (Klebsiella) must not match"
print("OK  neg-control E. coli flagged against the sample at {0:.1f}% identity".format(pid * 100))

# ---- scenario 2: neg control has its own unrelated organism -> CLEAN (no flag) -----------------
fasta2 = write_fasta(
    [
        ("barcode01|sample|0|Escherichia coli", ecoli),
        ("barcode03|sample|0|Klebsiella pneumoniae", kleb),
        ("barcode02|negative control|0|unclassified", staph),
    ]
)
recs2 = ccc.load_records(fasta2)
neg_hits2, _ = ccc.find_matches(recs2, min_identity=0.99)
assert neg_hits2 == [], "a neg control unrelated to every sample must not flag"
print("OK  neg control unrelated to all samples -> clean, no flag")

# ---- scenario 3: two DIFFERENT patients share E. coli -> sample<->sample carryover signal ------
fasta3 = write_fasta(
    [
        ("barcode01|sample|0|Escherichia coli", ecoli),
        ("barcode04|sample|0|Escherichia coli", mutate(ecoli, 2)),
        ("barcode05|sample|0|Klebsiella pneumoniae", kleb),
    ]
)
recs3 = ccc.load_records(fasta3)
neg_hits3, sample_hits3 = ccc.find_matches(recs3, min_identity=0.99)
assert neg_hits3 == [], "no negative control here"
assert len(sample_hits3) == 1, sample_hits3
print("OK  two samples sharing E. coli surfaced as a sample<->sample pair (possible carryover)")

# ---- end-to-end: main() writes the TSV and prints a summary -----------------------------------
import io, contextlib, argparse

buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    ccc.main(argparse.Namespace(fasta=fasta, min_identity=0.99, out=out))
printed = buf.getvalue()
assert "CONTAMINATION SIGNAL" in printed and "barcode02" in printed, printed
rows = open(out).read().splitlines()
assert rows[0].startswith("kind\t") and any(r.startswith("neg_vs_sample\t") for r in rows[1:]), rows
print("OK  main() writes a TSV and prints a human summary")

# ---- scenario 4: the --medaka-dir loader against a published-output tree -----------------------
root = tempfile.mkdtemp(prefix="medaka_pass_")


def write_consensus(barcode, cluster, s):
    d = os.path.join(root, "{0}_{1}_consensus_medaka".format(barcode, cluster))
    os.makedirs(d)
    with open(os.path.join(d, "consensus.fasta"), "w") as fh:
        fh.write(">{0}_{1}\n{2}\n".format(barcode, cluster, s))


write_consensus("barcode01", "0", ecoli)  # sample
write_consensus("barcode03", "0", kleb)  # sample
write_consensus("barcode02", "0", mutate(ecoli, 4))  # negative control, shares E. coli
recs4 = ccc.load_from_medaka_dir(root, ["barcode02"], [])
assert {(r["barcode"], r["status"]) for r in recs4} >= {
    ("barcode02", "negative control"),
    ("barcode01", "sample"),
}, recs4
neg_hits4, _ = ccc.find_matches(recs4, min_identity=0.99)
assert len(neg_hits4) == 1 and neg_hits4[0][0]["barcode"] == "barcode02", neg_hits4
print("OK  --medaka-dir loader labels controls from barcode and flags the match")

print("\nALL CROSS-CONTAMINATION CHECKS PASSED")
