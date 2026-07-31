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

1. Reads one combined FASTA whose headers encode the source: `>barcode|status|cluster|species`.
2. For every **negative-control** consensus, finds its most similar **sample** consensus by
   nucleotide identity (edlib infix alignment, so a partial consensus embedded in a longer one still
   matches).
3. Flags pairs at or above `--min-identity` (default **0.99**).
4. As a secondary signal, reports **sample↔sample** near-identical pairs from *different* barcodes
   (possible carryover / index hopping).
5. Writes `cross_contamination.tsv` and prints a human summary.

Tested on synthetic ~1500 bp sequences in `tests/cross_contamination_check.py`:
a neg-control consensus 99.7% identical to a sample is flagged; an unrelated neg control is not;
two patients sharing an organism surface as a sample↔sample pair.

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
- **Threshold is a policy choice.** 0.99 is a starting point; validate against real runs with known
  contamination before trusting it clinically.
- **Dependency:** the check needs `edlib` (tiny, pip-installable) in the report/clustering container,
  or swap the identity function for an external aligner (vsearch `--allpairs_global`, or minimap2).

## To validate / next steps

1. Run on real runs (including a known-contaminated one) and eyeball the flagged pairs + identities.
2. Decide the threshold and whether sample↔sample pairs should be reported.
3. Wire the TSV into the report (raise the neg-control light with the matched sample + identity).
4. Add `edlib` to the container (or switch to vsearch/minimap2) before enabling by default.
