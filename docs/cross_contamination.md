# Sequence-level cross-contamination check (prototype)

**Branch:** `contamination-clustering` — experimental, not wired into the default pipeline yet.

## Idea

Today's negative-control check compares **species names** between the negative control and a sample.
That is fragile: it depends on the classifier naming the same organism the same way in both, and on
abundance ranking. A cleaner signal is to compare the actual **consensus sequences**.

At the end of a run every barcode (samples, negative control, positive control) has produced one or
more polished consensus sequences. Pool them all and ask: *is any negative-control consensus
near-identical to a sample's consensus?* If so, the same DNA is in both — reagent/environmental
contamination or cross-contamination — no matter how each was named.

## How the prototype works

`bin/consensus_cross_contamination.py`:

1. Gathers all consensus, either from a combined FASTA (`--fasta`, headers `barcode|status|cluster|species`)
   or straight from a run's medaka output (`--medaka-dir <outdir>/medaka_pass`, with the control
   barcodes named via `--negative` / `--positive`).
2. Compares them **all-vs-all with minimap2** (`-c -X -x ava-ont`), which reports identity over the
   **aligned overlap** plus the overlap length.
3. Flags a **negative-control** consensus against its best **sample** match when identity ≥
   `--min-identity` (default **0.99**) **and** the overlap is ≥ `--min-overlap` bp (default **300**).
   Scoring identity over the overlap (not the whole sequence) means two consensus covering different,
   only-partially-overlapping stretches of 16S are compared fairly; the minimum-overlap rule stops a
   short, highly-conserved stretch from raising a false flag.
4. As a secondary signal, reports **sample↔sample** near-identical pairs from *different* barcodes
   (possible carryover / index hopping).
5. Writes `cross_contamination.tsv` (with an `overlap_bp` column) and prints a human summary.

If minimap2 isn't on `PATH`, it falls back to a built-in edit-distance identity (edlib, or stdlib
difflib) — install-free and fine for a quick look at full-length data, but it measures identity over
the whole shorter sequence, so it can miss offset/partial overlaps (which is exactly why minimap2 is
preferred).

Tested in `tests/cross_contamination_check.py`: the built-in path on synthetic ~1500 bp sequences,
and the minimap2 path via synthetic PAF (identity-over-overlap flags a match; the min-overlap rule
filters short high-identity hits; a 400 bp @ 99.5% offset overlap is caught).

## Proposed pipeline integration

A single collectFile in `workflows/nanopath.nf` builds the combined FASTA, then one process runs the
check (`modules/local/cross_contamination/main.nf`):

```groovy
// pool every barcode's consensus into one FASTA; headers encode barcode|status|cluster
ch_all_consensus = MEDAKA_PASS.out.consensus
    .map { meta, consensus, draft_log, cluster ->
        def seq = consensus.readLines().findAll { !it.startsWith('>') }.join('')
        ">${meta.id}|${meta.status}|${cluster}\n${seq}\n"
    }
    .collectFile(name: 'all_consensus.fasta', newLine: false)

CONSENSUS_CROSS_CHECK( ch_all_consensus )
```

(Species can be added to the header later by joining with the chosen-classifier call per cluster.)
The resulting `cross_contamination.tsv` can be published and/or surfaced in the run report, and a
flagged neg↔sample pair could raise the negative-control light with a specific, sequence-based reason.

## Caveats (important)

- **16S can't resolve strains.** Different strains of one species are often >99% identical, so a flag
  means "the same organism/sequence is in the negative control and a sample" — exactly what we want to
  surface — but it does not prove the *direction* of contamination. Report the identity: 100% is much
  stronger evidence of shared physical DNA than 99%.
- **Genuine shared species.** Two patients truly infected with the same organism will match at the
  sequence level; that is why the neg control (which should be sterile) is the anchor. Sample↔sample
  matches are a softer "worth a look", not a hard fail.
- **Thresholds are a policy choice.** identity 0.99 and overlap 300 bp are starting points; validate
  against real runs (with and without known contamination) and tune before trusting them clinically.
  The minimap2 `--preset` (default `ava-ont`) may also need checking for short 16S amplicons.
- **Dependency:** the preferred path needs `minimap2` on `PATH` (add it to the module's container).
  Without it the built-in edit-distance fallback runs, which does not handle offset/partial overlaps.

## To validate / next steps

1. Run on real runs (including a known-contaminated one) and eyeball the flagged pairs, identities,
   and overlap lengths.
2. Decide the identity/overlap thresholds and whether sample↔sample pairs should be reported.
3. Wire the TSV into the report (raise the neg-control light with the matched sample + identity).
4. Add `minimap2` to the container before enabling by default.
