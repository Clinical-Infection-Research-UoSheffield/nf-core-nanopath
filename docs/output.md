# nf-core/nanopath: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to the top-level results directory specified with `--outdir`.

## Pipeline overview

The pipeline processes Nanopore 16S/ITS2 amplicon reads through the following steps, each of which produces its own output directory:

1. _(optional)_ **CAT_FASTQS** — Concatenation of per-barcode FASTQ directories
2. **FASTP** — Read length filtering and quality control
3. **FastQC** — Read quality metrics
4. _(optional)_ **KRAKEN2** — Pre-filtering of unclassified reads
5. **KMER_FREQS** — 5-mer frequency computation
6. **READ_CLUSTERING** — UMAP dimensionality reduction + HDBSCAN clustering
7. **SPLIT_CLUSTERS** — Splitting reads by cluster assignment
8. **CANU_CORRECTION** — Error correction per cluster
9. **DRAFT_SELECTION** — Best representative read selection via FastANI
10. **RACON_PASS** — Racon polishing
11. **MEDAKA_PASS** — Medaka consensus polishing
12. **\*_CLASSIFICATION** — Taxonomic classification (Kraken2 / BLAST / SeqMatch / all)
13. **JOIN_RESULTS** — Joining per-cluster results into a per-sample table
14. **GET_ABUNDANCE** — Relative abundance calculation at multiple taxonomic levels ⭐
15. _(optional)_ **GENERATE_REPORTS** — Per-patient HTML clinical reports ⭐
16. **MultiQC** — Aggregated QC report
17. **Pipeline info** — Nextflow execution metadata

---

## Output directories

### `cat/`

> Only produced when `--fastq_dir` is provided.

<details markdown="1">
<summary>Output files</summary>

- `cat/`
  - `barcode*.fastq.gz`: One concatenated FASTQ file per barcode directory found under `--fastq_dir`.

</details>

All FASTQ files within each `barcode*` subdirectory are merged into a single file before any downstream processing.

---

### `fastp/`

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `*.fastq.gz`: Length-filtered reads (reads outside `--min_read_length`–`--max_read_length` are removed).
  - `*.fastp.json`: Per-sample JSON QC report (read counts before/after filtering, length distribution).
  - `*.fastp.html`: Per-sample HTML QC report.
  - `*.fastp.log`: FASTP log file.

</details>

[FASTP](https://github.com/OpenGene/fastp) is used solely for read length filtering (default 1400–1700 bp for 16S amplicons). Samples with zero reads remaining after filtering are automatically discontinued with a warning in the Nextflow log.

---

### `fastqc/`

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC quality report per sample.
  - `*_fastqc.zip`: Zip archive containing the report, tab-delimited data and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports are run on the **raw** input reads (before FASTP filtering) and are collected into the MultiQC report. They are useful for assessing overall read quality, length distribution and potential adapter contamination.

---

### `kmer_freqs/`

<details markdown="1">
<summary>Output files</summary>

- `kmer_freqs/`
  - `*_freqs.txt`: Tab-delimited table of log-normalised 5-mer frequencies per read. Columns: `read`, `length`, followed by one column per canonical 5-mer.

</details>

The k-mer frequency matrix is the input to the UMAP/HDBSCAN clustering step. Each row represents one read; values are log-transformed normalised counts of all canonical 5-mers (forward + reverse complement combined). Up to `--umap_set_size` reads are used (default 100,000).

---

### `read_clustering/`

<details markdown="1">
<summary>Output files</summary>

- `read_clustering/`
  - `*_hdbscan_output.tsv`: Tab-delimited cluster assignment table. Columns: `read`, `length`, `D1` (UMAP dimension 1), `D2` (UMAP dimension 2), `bin_id` (cluster ID; reads assigned to noise cluster −1 are excluded).
  - `*_hdbscan_output.png`: UMAP scatter plot coloured by cluster assignment.

</details>

UMAP reduces the k-mer frequency matrix to 2D, and HDBSCAN identifies clusters of reads. Each cluster typically corresponds to a distinct organism or strain present in the sample. The scatter plot provides a visual overview of the clustering result and the number of clusters detected.

Key parameters controlling this step: `--umap_n_neighbors`, `--umap_min_dist`, `--min_cluster_size`, `--min_samples`, `--cluster_sel_epsilon`.

---

### `canu_correction/`

<details markdown="1">
<summary>Output files</summary>

- `canu_correction/`
  - `*_cluster<N>.report`: Canu assembly report per cluster, summarising input read statistics and correction outcomes.

</details>

[Canu](https://github.com/marbl/canu) performs error correction on the reads in each cluster. The corrected reads (`.correctedReads.fasta.gz`) are passed to the next step but are not published by default. The `.report` file is useful for diagnosing clusters with low coverage or read count.

---

### `draft_selection/`

<details markdown="1">
<summary>Output files</summary>

- `draft_selection/`
  - `*_cluster<N>_draft_read.fasta`: Single FASTA record of the best representative (draft) read for each cluster, selected by highest mean FastANI identity to all other reads in the cluster.

</details>

[FastANI](https://github.com/ParBLiSS/FastANI) computes pairwise average nucleotide identity between all corrected reads in a cluster. The read with the highest mean ANI is selected as the draft template for polishing. Fragment length is set automatically based on assay type: 1200 bp for `16S`, 400 bp for `ITS2`.

---

### `racon_pass/`

<details markdown="1">
<summary>Output files</summary>

- `racon_pass/`
  - `*_cluster<N>_racon_consensus.fasta`: Racon-polished consensus sequence per cluster.

</details>

[Minimap2](https://github.com/lh3/minimap2) aligns all corrected reads to the draft read, and [Racon](https://github.com/isovic/racon) generates a polished consensus. If Racon fails due to insufficient overlaps, the draft read is used as the consensus and a warning is logged.

---

### `medaka_pass/`

<details markdown="1">
<summary>Output files</summary>

- `medaka_pass/`
  - `*_cluster<N>_consensus_medaka/consensus.fasta`: Final polished consensus sequence per cluster.

</details>

[Medaka](https://github.com/nanoporetech/medaka) further polishes the Racon consensus using a neural-network-based approach. If Medaka fails, the Racon consensus is copied as output. The `consensus.fasta` files from this step are the sequences submitted to taxonomic classification.

---

### `blast_classification/`, `seqmatch_classification/`, `kraken2_classification/`, `full_classification/`

The output directory name reflects the `--classification` parameter (`blast`, `seqmatch`, `kraken2`, or `full`).

<details markdown="1">
<summary>Output files</summary>

- `<classifier>_classification/`
  - `*_cluster<N>_*_consensus_classification.csv`: Top classification hits for the cluster consensus in semicolon-delimited format.
  - `*_cluster<N>_classification.log`: Log file carrying the cluster metadata and top classification result through to the `JOIN_RESULTS` step.
  - `*_cluster<N>_*_classification_out.tsv`: _(Kraken2 and full only)_ Full Kraken2 per-read output TSV.

</details>

**BLAST** (`--classification blast`): Runs `blastn` (megablast, no dust filter) against `--blast_db`. Output columns (semicolon-delimited): `scientific_name; taxid; evalue; length; percent_identity; bitscore`. Top 5 hits by identity/bitscore are retained.

**SeqMatch** (`--classification seqmatch`): Runs RDP `SequenceMatch` against `--seqmatch_db`, joined with `--seqmatch_accession` to recover taxonomy. Output columns: `scientific_name; taxid; seqmatch_score; name; species; genus; family; order`.

**Kraken2** (`--classification kraken2`): Runs `kraken2 --report` against `--kraken2_db`. If `--reclassifyOnFail` is set and the classification does not reach species level, SeqMatch is automatically used as a fallback.

**Full** (`--classification full`): All three classifiers run in parallel on the same consensus. The `get_abundance.py` script subsequently selects the best result per cluster based on classification depth (species-level calls preferred; otherwise the classifier with most non-null taxonomy fields is chosen).

---

### `join_results/`

<details markdown="1">
<summary>Output files</summary>

- `join_results/`
  - `<sample>.nanoclust_out.txt`: Semicolon-delimited text file with one row per cluster, joining metadata, classification results and (when `--taxonomy` is provided) full lineage annotation.

</details>

Column headers depend on `--classification`:

| Mode | Columns |
|------|---------|
| `kraken2` / `seqmatch` | `id; reads_in_cluster; used_for_consensus; reads_after_corr; draft_id; sciname; taxid; class_level; name; species; genus; family; order` |
| `blast` | `id; reads_in_cluster; used_for_consensus; reads_after_corr; draft_id; sciname; taxid; length; per_ident` |
| `full` | All columns from all three classifiers concatenated |

When `--taxonomy` is provided, the lineage fields (`name`, `species`, `genus`, `family`, `order`) are populated by a lookup against the taxonomy file using the reported taxid.

---

### `get_abundance/` ⭐ Primary result

> **These are the main quantitative outputs of the pipeline** — per-sample relative abundance tables at four taxonomic levels, directly suitable for downstream analysis and visualisation.

<details markdown="1">
<summary>Output files</summary>

- `get_abundance/`
  - `rel_abundance_<sample>_S.csv`: Relative abundance at **species** level. **This is the most commonly used output file.**
  - `rel_abundance_<sample>_G.csv`: Relative abundance at **genus** level.
  - `rel_abundance_<sample>_F.csv`: Relative abundance at **family** level.
  - `rel_abundance_<sample>_O.csv`: Relative abundance at **order** level.

</details>

Each CSV has three columns:

| Column | Description |
|--------|-------------|
| `taxid` | Taxonomic name (resolved from `--taxonomy` file or Unipept API fallback) |
| `rel_abundance` | Relative abundance as a percentage of total reads in the sample |
| `reads` | Absolute read count attributed to this taxon |

Results are sorted descending by relative abundance. For `--classification full`, the script automatically selects the best classifier result per cluster: species-level Kraken2 hits are preferred; otherwise the classifier with the most non-null taxonomy fields is chosen.

---

### `generate_reports/` ⭐ Primary result (clinical mode)

> **These are the primary end-user deliverable in clinical mode** — one self-contained HTML report per patient, ready to share without any additional tooling.

> Only produced when both `--clinical` and `--generateReports` are set.

<details markdown="1">
<summary>Output files</summary>

- `generate_reports/`
  - `patient_report_<specimen_number>.html`: Self-contained HTML report per patient sample, viewable in any web browser.

</details>

Reports are generated using a University of Sheffield branded template (`assets/UoS_report_template.html`) and include:

- **Sample metadata table**: specimen number, barcode, assay type, run ID and sequencing start time.
- **Results table**: top 3 detected species with relative abundance (%) and read counts.
- **Control summary**: results from positive and negative controls run in the same batch.
- **QC statistics**: total reads after quality control.
- **Discontinued samples**: samples with `Status = discontinued` in the samplesheet receive a dedicated "not detected" report page listing their metadata.

---

### `multiqc/`

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: A standalone HTML report aggregating QC metrics from all samples and pipeline steps.
  - `multiqc_data/`: Directory containing parsed statistics from each tool.
  - `multiqc_plots/`: Directory containing static images from the report.

</details>

[MultiQC](http://multiqc.info) collects FastQC results across all samples and summarises software versions used in the pipeline run. The report is useful for quickly identifying samples with quality issues.

---

### `pipeline_info/`

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_<timestamp>.html`: Nextflow execution report (per-process resource usage).
  - `execution_timeline_<timestamp>.html`: Gantt chart of process execution times.
  - `execution_trace_<timestamp>.txt`: Tab-delimited resource trace per task.
  - `pipeline_dag_<timestamp>.html`: Directed acyclic graph of the pipeline.
  - `software_versions.yml` / `*_versions.yml`: Collated software version YAML used for MultiQC methods description.
  - `samplesheet.valid.csv`: The validated and normalised samplesheet generated by `SAMPLESHEET_CHECK`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) generates execution reports automatically. These are invaluable for troubleshooting failed runs, reviewing resource utilisation and confirming which software versions were used.
