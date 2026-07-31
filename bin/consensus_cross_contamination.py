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

Sequences are compared all-vs-all with minimap2, which reports identity over the *aligned overlap*
plus the overlap length. A negative-control consensus is flagged against its best sample match when
identity >= --min-identity AND the overlap is at least --min-overlap bp. Measuring identity over the
overlap (not the whole sequence) means two consensus that cover different, only-partially-overlapping
stretches of 16S are still compared fairly; the minimum-overlap rule stops a short, highly-conserved
stretch (16S has very conserved regions) from raising a false flag. sample<->sample near-identical
pairs are reported as a secondary "possible carryover" signal.

Caveat worth remembering: 16S is highly conserved, so different strains of one species are often
>99% identical. A flag therefore means "the same organism/sequence is in the negative control and a
sample" -- exactly what we want to surface -- but it cannot by itself prove the direction of
contamination. Report the identity so a human can judge (100% is stronger than 99%).

    consensus_cross_contamination.py --medaka-dir <outdir>/medaka_pass --negative barcode02 \
        --min-identity 0.99 --min-overlap 300 --out cross_contamination.tsv
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile

# Preferred engine: minimap2 all-vs-all. It reports identity over the *aligned region* plus the
# alignment length, so two consensus that cover different, only-partially-overlapping stretches of
# 16S are compared over their overlap (not penalised for the non-overlapping ends), and we can
# require a minimum overlap length. If minimap2 isn't on PATH we fall back to a built-in edit-
# distance identity (edlib, or stdlib difflib) -- fast and install-free, but it measures identity
# over the whole shorter sequence, so it can miss offset/partial overlaps. See docs/cross_contamination.md.
try:
    import edlib

    _HAVE_EDLIB = True
except ImportError:
    import difflib

    _HAVE_EDLIB = False


def _have_minimap2():
    return shutil.which("minimap2") is not None


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


def load_from_medaka_dir(medaka_dir, negatives, positives):
    """Build records directly from a run's published medaka output -- no combined FASTA needed.

    Globs '<dir>/**/<barcode>_<cluster>_consensus_medaka/consensus.fasta' (as published to
    <outdir>/medaka_pass/). Barcode and cluster come from the directory name; status is set from
    the --negative / --positive barcode lists, everything else is a sample. This is the one-command
    way to test the check on an existing run.
    """
    negs = {b.strip() for b in negatives if b and b.strip()}
    poss = {b.strip() for b in positives if b and b.strip()}
    recs, seen = [], set()
    pattern = os.path.join(medaka_dir, "**", "*_consensus_medaka", "consensus.fasta")
    for path in sorted(glob.glob(pattern, recursive=True)):
        real = os.path.realpath(path)
        if real in seen:
            continue
        seen.add(real)
        dirname = os.path.basename(os.path.dirname(path))
        stem = dirname[: -len("_consensus_medaka")] if dirname.endswith("_consensus_medaka") else dirname
        barcode, _, cluster = stem.rpartition("_")
        barcode = barcode or stem
        status = NEG if barcode in negs else POS if barcode in poss else SAMPLE
        seq = "".join(s for _, s in parse_fasta(path))  # concatenate any records in the file
        if not seq:
            continue
        recs.append(
            {
                "barcode": barcode,
                "status": status,
                "cluster": cluster or "?",
                "species": "",
                "seq": seq.upper(),
                "name": "{0} {1} (cluster {2})".format(barcode, status, cluster or "?"),
            }
        )
    return recs


def identity_builtin(a, b):
    """Fallback nucleotide identity in [0,1] when minimap2 is unavailable.

    With edlib: infix ('HW') alignment of the shorter into the longer, so a partial consensus
    embedded in a longer one still matches over its length. Without edlib: the stdlib difflib ratio.
    Both measure identity over the whole shorter sequence, so offset/partial overlaps score low --
    that is the known limitation minimap2 fixes.
    """
    if not a or not b:
        return 0.0
    if _HAVE_EDLIB:
        q, t = (a, b) if len(a) <= len(b) else (b, a)
        res = edlib.align(q, t, mode="HW", task="distance")
        dist = res["editDistance"]
        if dist is None or dist < 0:
            return 0.0
        return max(0.0, 1.0 - dist / len(q))
    return difflib.SequenceMatcher(None, a, b, autojunk=False).ratio()


def _run_minimap2(recs, preset):
    """Run minimap2 all-vs-all on the pooled consensus; return PAF lines. Sequences are written with
    numeric names (index into recs) so the PAF maps straight back."""
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
    try:
        for i, r in enumerate(recs):
            tmp.write(">{0}\n{1}\n".format(i, r["seq"]))
        tmp.close()
        # -c base-level alignment (so identity is exact), -X skip self/dual mappings (all-vs-all)
        cmd = ["minimap2", "-c", "-X", "-x", preset, tmp.name, tmp.name]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if proc.returncode != 0:
            raise RuntimeError("minimap2 failed:\n" + proc.stderr[-800:])
        return proc.stdout.splitlines()
    finally:
        os.unlink(tmp.name)


def parse_paf(lines):
    """Yield (query_index, target_index, identity, overlap_len) from PAF lines.

    PAF col 10 = matching bases, col 11 = alignment block length (matches+mismatches+indels).
    identity = col10/col11 (identity over the aligned region); overlap length = col11.
    """
    for ln in lines:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        try:
            q, t = int(f[0]), int(f[5])
            matches, block = int(f[9]), int(f[10])
        except ValueError:
            continue
        if block <= 0:
            continue
        yield q, t, matches / block, block


def compute_pairs(recs, min_identity, min_overlap, paf_lines=None, preset="ava-ont"):
    """Best (identity, overlap) for each cross-barcode pair meeting BOTH thresholds.

    Uses minimap2 (given paf_lines, or by running it when on PATH); otherwise the built-in
    edit-distance fallback. Returns {(i, j): (identity, overlap)} with i < j.
    """
    best = {}

    def consider(i, j, ident, ov):
        if i == j or recs[i]["barcode"] == recs[j]["barcode"]:
            return  # self, or same barcode (different clusters of one sample) -> not contamination
        if ident < min_identity or ov < min_overlap:
            return
        key = (min(i, j), max(i, j))
        if key not in best or ident > best[key][0]:
            best[key] = (ident, ov)

    if paf_lines is not None or _have_minimap2():
        if paf_lines is None:
            paf_lines = _run_minimap2(recs, preset)
        for q, t, ident, ov in parse_paf(paf_lines):
            consider(q, t, ident, ov)
    else:
        for i in range(len(recs)):
            for j in range(i + 1, len(recs)):
                ov = min(len(recs[i]["seq"]), len(recs[j]["seq"]))
                consider(i, j, identity_builtin(recs[i]["seq"], recs[j]["seq"]), ov)
    return best


def find_matches(recs, min_identity, min_overlap=300, paf_lines=None, preset="ava-ont"):
    """Flagged pairs as (rec_a, rec_b, identity, overlap): every negative-control consensus vs its
    best sample match, plus sample<->sample near-identical pairs. Thresholds: identity AND overlap."""
    best = compute_pairs(recs, min_identity, min_overlap, paf_lines, preset)
    neg_best = {}  # neg index -> (sample_rec, identity, overlap)
    sample_hits = []
    for (i, j), (ident, ov) in best.items():
        a, b = recs[i], recs[j]
        if a["status"] == NEG and b["status"] == SAMPLE:
            ni, srec = i, b
        elif b["status"] == NEG and a["status"] == SAMPLE:
            ni, srec = j, a
        elif a["status"] == SAMPLE and b["status"] == SAMPLE:
            sample_hits.append((a, b, ident, ov))
            continue
        else:
            continue
        if ni not in neg_best or ident > neg_best[ni][1]:
            neg_best[ni] = (srec, ident, ov)
    neg_hits = [(recs[ni], srec, ident, ov) for ni, (srec, ident, ov) in neg_best.items()]
    return neg_hits, sample_hits


def write_tsv(path, neg_hits, sample_hits):
    with open(path, "w") as fh:
        fh.write("kind\tidentity\toverlap_bp\tsource_a\tspecies_a\tsource_b\tspecies_b\n")
        for kind, hits in (("neg_vs_sample", neg_hits), ("sample_vs_sample", sample_hits)):
            for a, b, pid, ov in hits:
                fh.write(
                    "{0}\t{1:.4f}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(
                        kind, pid, ov, a["name"], a["species"], b["name"], b["species"]
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
        for a, b, pid, ov in sorted(neg_hits, key=lambda x: -x[2]):
            lines.append(
                "  - {0} is {1:.1f}% identical to {2} (over {3} bp)".format(a["name"], pid * 100, b["name"], ov)
            )
    elif n_neg:
        lines.append("No negative-control consensus matched any sample above threshold. Clean.")
    if sample_hits:
        lines.append("Also: {0} near-identical sample<->sample pair(s) (possible carryover):".format(len(sample_hits)))
        for a, b, pid, ov in sorted(sample_hits, key=lambda x: -x[2]):
            lines.append("  - {0} <-> {1} : {2:.1f}% (over {3} bp)".format(a["name"], b["name"], pid * 100, ov))
    return "\n".join(lines)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--fasta", help="combined FASTA of all consensus, headers 'barcode|status|cluster|species'")
    src.add_argument(
        "--medaka-dir",
        help="a run's medaka output dir (e.g. <outdir>/medaka_pass); consensus are found and labelled automatically",
    )
    p.add_argument(
        "--negative", default="", help="comma-separated barcode(s) that are negative controls (with --medaka-dir)"
    )
    p.add_argument(
        "--positive", default="", help="comma-separated barcode(s) that are positive controls (with --medaka-dir)"
    )
    p.add_argument("--min-identity", type=float, default=0.99, help="identity over the overlap to flag a match [0.99]")
    p.add_argument("--min-overlap", type=int, default=300, help="minimum aligned overlap in bp to flag a match [300]")
    p.add_argument("--preset", default="ava-ont", help="minimap2 preset for the all-vs-all alignment [ava-ont]")
    p.add_argument("--out", default="cross_contamination.tsv", help="output TSV of flagged pairs")
    return p.parse_args()


def main(args):
    if _have_minimap2():
        sys.stderr.write(
            "Using minimap2 (identity over the aligned overlap; min overlap {0} bp).\n".format(args.min_overlap)
        )
    else:
        sys.stderr.write(
            "NOTE: minimap2 not found; using the built-in {0} fallback, which measures identity over "
            "the whole shorter sequence and may miss offset/partial overlaps.\n".format(
                "edlib" if _HAVE_EDLIB else "difflib"
            )
        )
    if getattr(args, "medaka_dir", None):
        recs = load_from_medaka_dir(args.medaka_dir, args.negative.split(","), args.positive.split(","))
    else:
        recs = load_records(args.fasta)
    n_neg = sum(1 for r in recs if r["status"] == NEG)
    n_sample = sum(1 for r in recs if r["status"] == SAMPLE)
    neg_hits, sample_hits = find_matches(recs, args.min_identity, args.min_overlap, preset=args.preset)
    write_tsv(args.out, neg_hits, sample_hits)
    print(summarise(neg_hits, sample_hits, n_neg, n_sample))
    return neg_hits


if __name__ == "__main__":
    main(parse_args())
