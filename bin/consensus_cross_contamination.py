#!/usr/bin/env python
"""Cross-contamination check by pooling all consensus sequences of a run.

PROTOTYPE (experiment on the contamination-clustering branch).

Rather than matching *species names* between the negative control and a sample (fragile: depends on
taxonomy naming and abundance ranking), this compares the actual consensus *sequences*. If a
negative-control consensus is near-identical to a sample's consensus, that is direct evidence of
shared DNA -- reagent/environmental contamination or cross-contamination -- however the classifier
happened to name them.

Input: one combined FASTA whose headers encode the source, pipe-separated:
    >{barcode}|{status}|{cluster}|{species}
       barcode : e.g. barcode01
       status  : 'sample' | 'negative control' | 'positive control'
       cluster : cluster id within that barcode
       species : optional, for nicer reporting (may be empty/omitted)

For every negative-control sequence we find its most similar *sample* sequence (nucleotide identity
via edlib infix alignment) and flag pairs at or above --min-identity. We also report sample<->sample
near-identical pairs (possible cross-contamination / index hopping) as a secondary signal.

Caveat worth remembering: 16S is highly conserved, so different strains of one species are often
>99% identical. A flag therefore means "the same organism/sequence is in the negative control and a
sample" -- exactly what we want to surface -- but it cannot by itself prove the direction of
contamination. Report the identity so a human can judge (100% is stronger than 99%).

    consensus_cross_contamination.py --fasta all_consensus.fasta --min-identity 0.99 \
        --out cross_contamination.tsv
"""

import argparse
import sys

try:
    import edlib
except ImportError:  # pragma: no cover - environment guard
    sys.stderr.write("ERROR: this check needs the 'edlib' package (pip install edlib).\n")
    raise

SAMPLE = "sample"
NEG = "negative control"
POS = "positive control"


def parse_fasta(path):
    """Yield (header, sequence) pairs from a FASTA file (uppercased, newlines stripped)."""
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            elif line:
                seq.append(line.strip().upper())
    if header is not None:
        yield header, "".join(seq)


def parse_label(header):
    """Split a '>barcode|status|cluster|species' header into a record dict."""
    parts = (header.split("|") + ["", "", "", ""])[:4]
    barcode, status, cluster, species = (p.strip() for p in parts)
    return {
        "barcode": barcode or "?",
        "status": (status or SAMPLE).lower(),
        "cluster": cluster or "?",
        "species": species,
    }


def load_records(path):
    recs = []
    for header, seq in parse_fasta(path):
        if not seq:
            continue
        r = parse_label(header)
        r["seq"] = seq
        r["name"] = "{0} {1} (cluster {2})".format(r["barcode"], r["status"], r["cluster"])
        recs.append(r)
    return recs


def identity(a, b):
    """Nucleotide identity in [0,1] via edlib infix ('HW') alignment of the shorter into the longer.

    Infix alignment tolerates length differences (a partial consensus embedded in a longer one still
    matches over its length), which is the realistic contamination case.
    """
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    if not q:
        return 0.0
    res = edlib.align(q, t, mode="HW", task="distance")
    dist = res["editDistance"]
    if dist is None or dist < 0:
        return 0.0
    return max(0.0, 1.0 - dist / len(q))


def find_matches(recs, min_identity):
    """Return flagged (a, b, identity) pairs: every negative-control seq vs its best sample seq,
    plus sample<->sample near-identical pairs. Only pairs at/above min_identity are returned."""
    neg = [r for r in recs if r["status"] == NEG]
    sample = [r for r in recs if r["status"] == SAMPLE]

    neg_hits = []
    for n in neg:
        best, best_id = None, 0.0
        for s in sample:
            pid = identity(n["seq"], s["seq"])
            if pid > best_id:
                best, best_id = s, pid
        if best is not None and best_id >= min_identity:
            neg_hits.append((n, best, best_id))

    sample_hits = []
    for i in range(len(sample)):
        for j in range(i + 1, len(sample)):
            if sample[i]["barcode"] == sample[j]["barcode"]:
                continue  # same barcode, different cluster -> expected, not cross-contamination
            pid = identity(sample[i]["seq"], sample[j]["seq"])
            if pid >= min_identity:
                sample_hits.append((sample[i], sample[j], pid))

    return neg_hits, sample_hits


def write_tsv(path, neg_hits, sample_hits):
    with open(path, "w") as fh:
        fh.write("kind\tidentity\tsource_a\tspecies_a\tsource_b\tspecies_b\n")
        for a, b, pid in neg_hits:
            fh.write(
                "neg_vs_sample\t{0:.4f}\t{1}\t{2}\t{3}\t{4}\n".format(
                    pid, a["name"], a["species"], b["name"], b["species"]
                )
            )
        for a, b, pid in sample_hits:
            fh.write(
                "sample_vs_sample\t{0:.4f}\t{1}\t{2}\t{3}\t{4}\n".format(
                    pid, a["name"], a["species"], b["name"], b["species"]
                )
            )


def summarise(neg_hits, sample_hits, n_neg, n_sample):
    lines = []
    if n_neg == 0:
        lines.append("No negative-control consensus sequences to compare.")
    if neg_hits:
        lines.append(
            "CONTAMINATION SIGNAL: {0} negative-control consensus/consensuses match a sample:".format(len(neg_hits))
        )
        for a, b, pid in sorted(neg_hits, key=lambda x: -x[2]):
            lines.append("  - {0} is {1:.1f}% identical to {2}".format(a["name"], pid * 100, b["name"]))
    elif n_neg:
        lines.append("No negative-control consensus matched any sample above threshold. Clean.")
    if sample_hits:
        lines.append("Also: {0} near-identical sample<->sample pair(s) (possible carryover):".format(len(sample_hits)))
        for a, b, pid in sorted(sample_hits, key=lambda x: -x[2]):
            lines.append("  - {0} <-> {1} : {2:.1f}%".format(a["name"], b["name"], pid * 100))
    return "\n".join(lines)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--fasta", required=True, help="combined FASTA of all consensus, headers 'barcode|status|cluster|species'"
    )
    p.add_argument("--min-identity", type=float, default=0.99, help="identity threshold to flag a match [0.99]")
    p.add_argument("--out", default="cross_contamination.tsv", help="output TSV of flagged pairs")
    return p.parse_args()


def main(args):
    recs = load_records(args.fasta)
    n_neg = sum(1 for r in recs if r["status"] == NEG)
    n_sample = sum(1 for r in recs if r["status"] == SAMPLE)
    neg_hits, sample_hits = find_matches(recs, args.min_identity)
    write_tsv(args.out, neg_hits, sample_hits)
    print(summarise(neg_hits, sample_hits, n_neg, n_sample))
    return neg_hits


if __name__ == "__main__":
    main(parse_args())
