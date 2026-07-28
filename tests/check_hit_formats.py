#!/usr/bin/env python
"""Level-1 check: run the real hit-details parsers on REAL classifier output.

Level 0 (tests/smoke_hit_details.py) proves the parsing logic is correct *given*
that the fixture formats are faithful. This script closes that gap: point it at a
directory containing genuine *_consensus_classification.csv files (a Nextflow work
dir, or a published results dir) and it runs the actual parsers from
bin/results_report.py over them and sanity-checks the extracted numbers.

    ./venv/bin/python tests/check_hit_formats.py <dir>        # searches recursively

Reports one line per file, then a summary. Exits non-zero if any file fails a hard
check (a strong signal the column layout no longer matches the parser). "WARN" lines
are things to eyeball, not necessarily bugs.
"""

import sys, os, types, glob, importlib.util

# stub heavy imports and load the real parsers (same trick as Level 0)
for name in ("aplanat", "aplanat.report", "bokeh", "bokeh.resources", "requests"):
    sys.modules[name] = types.ModuleType(name)
sys.modules["bokeh.resources"].INLINE = None
sys.modules["aplanat"].report = sys.modules["aplanat.report"]

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO_ROOT, "bin")
spec = importlib.util.spec_from_file_location("results_report", os.path.join(BIN, "results_report.py"))
rr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rr)


def _num(x):
    try:
        return float(str(x).strip())
    except (ValueError, TypeError):
        return None


def check_file(path):
    """Return (errors, warns) for a single per-cluster classification CSV."""
    classifier = rr.detect_classifier(path)
    errors, warns = [], []
    if classifier is None:
        return ["unrecognised filename (no _blastn_/_seqmatch_/_kraken2_ token)"], []

    if classifier == "blastn":
        recs = rr.blast_records(path)
        real = [r for r in recs if not rr._is_unclassified(r["species"])]
        for r in real:
            pid, bit = _num(r["pct_identity"]), _num(r["bitscore"])
            if pid is None or not (0 <= pid <= 100):
                errors.append(
                    f"pct_identity not a 0-100 number ({r['pct_identity']!r}) -> column shift? "
                    f"(species={r['species']!r}) -- likely a ';' inside sscinames"
                )
            if bit is None:
                warns.append(f"bitscore not numeric ({r['bitscore']!r})")
        print(f"  BLAST    {os.path.basename(path):55} hits={len(recs)} top={recs[0]['species'] if recs else '-'}")

    elif classifier == "seqmatch":
        recs = rr.seqmatch_records(path)
        real = [r for r in recs if not rr._is_unclassified(r["species"])]
        scores = [_num(r["s_ab_score"]) for r in real]
        for r, s in zip(real, scores):
            if s is None or not (0 <= s <= 1.0001):
                errors.append(
                    f"s_ab_score not a 0-1 number ({r['s_ab_score']!r}) -> column shift? " f"(species={r['species']!r})"
                )
        if scores and any(s is not None for s in scores):
            ok = [s for s in scores if s is not None]
            if ok != sorted(ok, reverse=True):
                warns.append("s_ab scores not in descending order (expected best-first)")
        print(f"  SeqMatch {os.path.basename(path):55} hits={len(recs)} top={recs[0]['species'] if recs else '-'}")

    elif classifier == "kraken2":
        recs = rr.kraken2_records(path)
        r = recs[0]
        if not str(r["species"]).strip():
            errors.append("kraken2 assigned name is empty -> report layout not tab-separated 6-col?")
        pct = _num(r["lca_reads_pct"])
        if r["lca_reads_pct"] != "" and (pct is None or not (0 <= pct <= 100)):
            errors.append(f"reads %% not a 0-100 number ({r['lca_reads_pct']!r})")
        print(f"  Kraken2  {os.path.basename(path):55} lca={r['species']!r} reads%={r['lca_reads_pct']}")

    return errors, warns


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("ERROR: give a directory to search for *_consensus_classification.csv")
        return 2
    root = sys.argv[1]
    files = sorted(glob.glob(os.path.join(root, "**", "*_consensus_classification.csv"), recursive=True))
    if not files:
        print(f"No *_consensus_classification.csv found under {root}")
        print("Tip: point at a Nextflow work/ dir or the published {blast,seqmatch,kraken2}_classification/ dirs.")
        return 2

    print(f"checking {len(files)} file(s) under {root}\n")
    total_err, total_warn = 0, 0
    for f in files:
        try:
            errors, warns = check_file(f)
        except Exception as e:  # a parse blow-up is itself a format-mismatch signal
            errors, warns = [f"parser raised {type(e).__name__}: {e}"], []
        for e in errors:
            print(f"    FAIL {os.path.basename(f)}: {e}")
            total_err += 1
        for wtext in warns:
            print(f"    WARN {os.path.basename(f)}: {wtext}")
            total_warn += 1

    print(f"\n{'FAILED' if total_err else 'OK'}: {len(files)} files, {total_err} error(s), {total_warn} warning(s)")
    return 1 if total_err else 0


if __name__ == "__main__":
    sys.exit(main())
