# ![nf-core/nanopath](docs/images/nf-core-nanopath_logo_light.png#gh-light-mode-only) ![nf-core/nanopath](docs/images/nf-core-nanopath_logo_dark.png#gh-dark-mode-only)

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/nanopath/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**nf-core/nanopath** pipeline is a bioinformatics tool designed for the processing of nanopore 16S/ITS sequencing data. It employs advanced clustering methods to group similar genetic sequences, facilitating the classification and reporting of bacterial and fungal constituents within input samples. The utilization of these clustering techniques contributes to a reduction in noise, minimizing the occurrence of false negatives and positives in the results. The pipeline thus serves as a reliable resource for obtaining precise insights into the microbial composition of the analyzed samples especially wihin the clinical setting. 

# ![NanopathPipeline](docs/images/16S_Pipeline.png#gh-light-mode-only) ![NanopathPipeline](docs/images/16S_Pipeline_darkmode.png#gh-dark-mode-only)

1. **Initialize the data:**
      - If a fastq directory is provided:
         - Concatenate fastq files using CAT_FASTQS.

2. **Validate input:**
      - Use the INPUT_CHECK subworkflow to read samplesheet, validate, and stage input files.
      - Branch reads based on their status (discontinued or samples).

3. **Perform Quality Control:**
      - Run ([`FASTP`](https://github.com/OpenGene/fastp)) for quality control, filtering, and preprocessing.
      - Filter out samples with no reads left after FASTP.
      - Run ([`FASTQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) on the processed reads.

4. **Classify and Cluster:**
      - If specified, remove unclassified reads using ([`KRAKEN2`](https://github.com/DerrickWood/kraken2)).
      - Subset reads based on specified parameters (default 100k reads to keep memory requirements reasonable).
      - Perform k-mer frequency analysis with KMER_FREQS.
      - Perform read clustering with READ_CLUSTERING using ([`HDBSCAN`](https://github.com/scikit-learn-contrib/hdbscan)) and ([`UMAP`](https://umap-learn.readthedocs.io/en/latest/)).

5. **Split Clusters and Correct Errors:**
      - Split clusters.
      - Perform error correction using ([`CANU`](https://github.com/marbl/canu)).

6. **Select and Polish Draft:**
      - Select draft reads using ([`FASTANI`](https://github.com/ParBLiSS/FastANI)).
      - Polish drafts using ([`RACON`](https://github.com/isovic/racon)).
      - Generate final consensus using ([`MEDAKA`](https://github.com/nanoporetech/medaka)).

7. **Classify Taxonomically:**
      - Based on chosen tool, classify consensus sequences with ([`BLAST`](https://www.ncbi.nlm.nih.gov/books/NBK279690/)), ([`SEQMATCH`](https://github.com/rdpstaff/SequenceMatch)), ([`KRAKEN`](https://github.com/DerrickWood/kraken2)) or all of them. 
      - Join classification results using JOIN_RESULTS.

8. **Estimate Abundace:**
      - Estimate abundance per sample per detected species. 

9. **Generate Reports:**
      - If report generation is chosen:
         - Generate HTML reports.

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

### Samplesheet preparation

The pipeline accepts three input modes. Choose the one that matches your data setup.

#### Standard mode

Prepare a CSV samplesheet with a `sample` column and a `fastq_1` column pointing to single-end Nanopore reads (`.fastq`, `.fastq.gz`, `.fq`, or `.fq.gz`):

```csv
sample,fastq_1
SAMPLE1,/path/to/SAMPLE1.fastq.gz
SAMPLE2,/path/to/SAMPLE2.fastq.gz
```

Run with:

```bash
nextflow run nf-core/nanopath \
   -profile <docker/singularity/conda> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --classification kraken2 \
   --kraken2_db /path/to/kraken2_db \
   --taxonomy /path/to/taxonomy.tsv
```

#### FastQ directory mode

If your reads are stored in per-barcode subdirectories (e.g. directly from a GridION/MinION run before manual concatenation), supply `--fastq_dir` pointing to the parent directory and a samplesheet with `sample` and `filename` columns:

```csv
sample,filename
SAMPLE1,barcode01
SAMPLE2,barcode02
```

The `CAT_FASTQS` step will automatically find and concatenate all FASTQ files inside each `barcode*` subdirectory.

```bash
nextflow run nf-core/nanopath \
   -profile <docker/singularity/conda> \
   --input samplesheet.csv \
   --fastq_dir /path/to/run_dir/fastq_pass \
   --outdir <OUTDIR> \
   --classification kraken2 \
   --kraken2_db /path/to/kraken2_db \
   --taxonomy /path/to/taxonomy.tsv
```

#### Clinical mode

For clinical sequencing runs with barcoded patient samples, use `--clinical` with a samplesheet containing `Specimen Number`, `Barcode`, `Status`, and `Assay` columns:

```csv
Specimen Number,Barcode,Status,Assay
PATIENT001,barcode01,active,16S
PATIENT002,barcode02,active,ITS2
PATIENT003,barcode03,discontinued,16S
```

- **Barcode** must follow the format `barcode01`–`barcode100`.
- **Status = discontinued** samples receive a report page noting no result, but classification is skipped. Any other status value is treated as normal.
- **Assay** must be `16S` (bacterial) or `ITS2` (fungal); this drives fragment length settings in the consensus polishing step.

Clinical mode is typically combined with `--fastq_dir` (or `--onGridion`) and `--generateReports`:

```bash
nextflow run nf-core/nanopath \
   -profile <docker/singularity/conda> \
   --input clinical_samplesheet.csv \
   --fastq_dir /path/to/run_dir/fastq_pass \
   --clinical \
   --generateReports \
   --outdir <OUTDIR> \
   --classification kraken2 \
   --kraken2_db /path/to/kraken2_db \
   --taxonomy /path/to/taxonomy.tsv
```

### Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--classification` | `kraken2` | Classification tool: `kraken2`, `blast`, `seqmatch`, or `full` (all three) |
| `--kraken2_db` | — | Path to a Kraken2 database directory |
| `--blast_db` | — | Path to a BLAST database prefix (e.g. `/path/to/db/16S_ribosomal_RNA`) |
| `--seqmatch_db` | — | Path to an RDP SeqMatch database file |
| `--seqmatch_accession` | — | Path to an RDP accession-to-taxonomy mapping file |
| `--taxonomy` | — | Path to a tab-delimited taxonomy file used to annotate classification results |
| `--remove_unclassified` | `false` | Pre-filter reads with Kraken2 before clustering to remove non-target sequences |
| `--reclassifyOnFail` | `false` | Fall back to SeqMatch if Kraken2 does not reach species level |
| `--umap_set_size` | `100000` | Max reads subsetted per sample for UMAP/HDBSCAN clustering |
| `--min_read_length` | `1400` | Minimum read length filter applied by FASTP (bp) |
| `--max_read_length` | `1700` | Maximum read length filter applied by FASTP (bp) |
| `--avg_amplicon_size` | `1.5k` | Expected amplicon size passed to Canu for error correction |
| `--clinical` | `false` | Enable clinical mode (alternative samplesheet format, per-patient reports) |
| `--generateReports` | `false` | Generate per-sample HTML patient reports (requires `--clinical`) |
| `--onGridion` | `false` | Auto-discover run metadata (`kit`, `run_id`, `seq_start`) from GridION output files located one level above `--fastq_dir` |
| `--report_file` | — | Manual alternative to `--onGridion`: path to the GridION `report*.html` file. Must be provided together with `--summary_file`. |
| `--summary_file` | — | Manual alternative to `--onGridion`: path to the GridION `final_summary*.txt` file. Must be provided together with `--report_file`. Both files are parsed to extract the sequencing kit, run ID and start time used in patient reports. |

For advanced clustering, polishing and Canu parameters, see the [full usage documentation](docs/usage.md).

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Pipeline output

All results are written to the directory specified with `--outdir`. The table below summarises the main output directories. For full file-level documentation see [docs/output.md](docs/output.md).

| Directory | Description |
|-----------|-------------|
| `cat/` | Per-barcode concatenated FASTQ files (only produced when `--fastq_dir` is used) |
| `fastqc/` | FastQC HTML quality reports and zip archives on the raw input reads |
| `fastp/` | Length-filtered reads (`.fastq.gz`) and per-sample JSON QC summaries |
| `read_clustering/` | HDBSCAN cluster assignment tables (`*_hdbscan_output.tsv`) and UMAP scatter plots (`*_hdbscan_output.png`) |
| `canu_correction/` | Canu error-correction reports per cluster |
| `draft_selection/` | Best representative (draft) read FASTA per cluster, selected by FastANI |
| `racon_pass/` | Racon-polished consensus FASTA per cluster |
| `medaka_pass/` | Final Medaka-polished consensus FASTA per cluster |
| `*_classification/` | Per-cluster classification output CSV and log files (directory name reflects chosen classifier); useful for detailed inspection of individual cluster results |
| `join_results/` | Per-sample joined classification table (`*.nanoclust_out.txt`) combining all cluster results |
| **`get_abundance/`** ⭐ | **Primary quantitative output** — relative abundance tables at species (`*_S.csv`), genus (`*_G.csv`), family (`*_F.csv`) and order (`*_O.csv`) level |
| **`generate_reports/`** ⭐ | **Per-patient clinical HTML reports** (`patient_report_*.html`; clinical mode only — requires `--clinical --generateReports`) |
| `multiqc/` | Aggregated MultiQC HTML report covering FastQC results and pipeline software versions |
| `pipeline_info/` | Nextflow execution report, timeline, trace, DAG and collated software version YAML |

## Credits

nf-core/nanopath was originally written by Magdalena Dabrowska as a project within the Universiy of Sheffield.

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/nanopath for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->
> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
