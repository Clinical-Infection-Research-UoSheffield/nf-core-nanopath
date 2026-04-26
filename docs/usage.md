# nf-core/nanopath: Usage

## Introduction

nf-core/nanopath processes single-end Nanopore 16S rRNA or ITS2 amplicon sequencing data. It filters, clusters, and polishes reads to build per-species consensus sequences, which are then taxonomically classified and summarised as relative abundance tables. An optional clinical reporting mode generates per-patient HTML reports.

## Input modes

The pipeline supports three mutually compatible input modes.

---

### 1. Standard mode

The simplest mode: provide a CSV samplesheet where each row points to one (single-end) FASTQ file per sample.

```csv
sample,fastq_1
SAMPLE1,/absolute/path/to/SAMPLE1.fastq.gz
SAMPLE2,/absolute/path/to/SAMPLE2.fastq.gz
```

| Column | Description |
|--------|-------------|
| `sample` | Unique sample identifier. Spaces are converted to underscores. |
| `fastq_1` | Absolute path to a `.fastq`, `.fastq.gz`, `.fq`, or `.fq.gz` file. |

Run command:

```bash
nextflow run nf-core/nanopath \
   -profile docker \
   --input samplesheet.csv \
   --outdir results/ \
   --classification kraken2 \
   --kraken2_db /path/to/kraken2_db \
   --taxonomy /path/to/taxonomy.tsv
```

---

### 2. FastQ directory mode

Use this when reads are organised in per-barcode subdirectories as produced by a GridION or MinION run (e.g. `fastq_pass/barcode01/*.fastq.gz`). The pipeline concatenates all files within each subdirectory before processing.

Samplesheet columns: `sample` and `filename` (the barcode directory name).

```csv
sample,filename
SAMPLE1,barcode01
SAMPLE2,barcode02
```

Pass `--fastq_dir` pointing to the parent directory that contains the `barcode*` subdirectories:

```bash
nextflow run nf-core/nanopath \
   -profile docker \
   --input samplesheet.csv \
   --fastq_dir /path/to/run_dir/fastq_pass \
   --outdir results/ \
   --classification kraken2 \
   --kraken2_db /path/to/kraken2_db \
   --taxonomy /path/to/taxonomy.tsv
```

Concatenated files are written to `<outdir>/cat/` before downstream processing.

---

### 3. Clinical mode

Designed for clinical diagnostic workflows with barcoded patient samples. Enable with `--clinical`.

Samplesheet columns: `Specimen Number`, `Barcode`, `Status`, `Assay`, plus any additional metadata columns you wish to propagate to patient reports.

```csv
Specimen Number,Barcode,Status,Assay
PATIENT001,barcode01,active,16S
PATIENT002,barcode02,active,ITS2
PATIENT003,barcode03,discontinued,16S
```

| Column | Description |
|--------|-------------|
| `Specimen Number` | Clinical specimen identifier (propagated to patient reports). |
| `Barcode` | Sequencing barcode in `barcode01`–`barcode100` format (case-insensitive). |
| `Status` | Any value other than `discontinued` will be processed normally. Use `discontinued` to skip classification for a sample and generate a "not detected" report page instead. |
| `Assay` | `16S` (bacterial) or `ITS2` (fungal). Drives fragment length in consensus polishing. |

Clinical mode is typically combined with `--fastq_dir` or `--onGridion`, and `--generateReports`:

```bash
nextflow run nf-core/nanopath \
   -profile docker \
   --input clinical_samplesheet.csv \
   --fastq_dir /path/to/run_dir/fastq_pass \
   --clinical \
   --generateReports \
   --outdir results/ \
   --classification kraken2 \
   --kraken2_db /path/to/kraken2_db \
   --taxonomy /path/to/taxonomy.tsv
```

#### GridION metadata

Run metadata (sequencing kit, run ID, start time) can be sourced in three ways:

| Method | Parameters |
|--------|------------|
| Auto-discovery (GridION) | `--onGridion` — automatically locates `report*.html` and `final_summary*.txt` one directory level above `--fastq_dir` |
| Manual file paths | `--report_file <path/to/report*.html>` **and** `--summary_file <path/to/final_summary*.txt>` — use this when the files are not in the expected location relative to `--fastq_dir`. Both must be provided together; they are parsed to extract the sequencing kit, run ID and sequencing start time for patient reports. |
| Manual values | `--kit <string>` `--run_id <string>` `--seq_start <string>` — use this to supply the values directly without providing any files. |

---

## Parameters

### Input / output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | **required** | Path to input samplesheet (CSV). |
| `--outdir` | **required** | Directory for all output files. |
| `--fastq_dir` | `null` | Directory containing per-barcode FASTQ subdirectories. |
| `--clinical` | `false` | Enable clinical samplesheet format and clinical report generation. |
| `--onGridion` | `false` | Auto-discover run metadata from GridION output files. |
| `--report_file` | `null` | Path to GridION HTML report file. |
| `--summary_file` | `null` | Path to GridION `final_summary*.txt` file. |
| `--kit` | `unknown` | Sequencing / barcoding kit name (used in patient reports). |
| `--run_id` | `unknown` | Run identifier (used in patient reports). |
| `--seq_start` | `unknown` | Sequencing start time string (used in patient reports). |
| `--generateReports` | `false` | Generate per-sample HTML patient reports (requires `--clinical`). |

### Quality control and read filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_read_length` | `1400` | Minimum read length (bp) passed to FASTP (`-l`). |
| `--max_read_length` | `1700` | Maximum read length (bp) passed to FASTP (`--length_limit`). |

Reads that fail FASTP length filtering are discarded. Samples with zero reads remaining after filtering are automatically discontinued with a warning.

### UMAP / HDBSCAN clustering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--umap_set_size` | `100000` | Maximum number of reads per sample used for UMAP dimensionality reduction. Reads beyond this limit are subsetted (first N reads). |
| `--umap_n_neighbors` | `15` | UMAP `n_neighbors` — size of local neighbourhood used to learn the manifold. |
| `--umap_min_dist` | `0.1` | UMAP `min_dist` — minimum distance between points in the 2D embedding. |
| `--min_cluster_size` | `50` | HDBSCAN minimum number of reads to form an independent cluster. |
| `--min_samples` | `n` | HDBSCAN `min_samples` — measure of clustering conservatism. `n` sets it to `None` (default behaviour). |
| `--cluster_sel_epsilon` | `0.5` | HDBSCAN `cluster_selection_epsilon` — minimum distance to separate clusters. |
| `--throughput` | `standard` | Compute resource label for clustering: `standard`, `high` (more memory), or `low`. |

### Consensus polishing

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--avg_amplicon_size` | `1.5k` | Expected amplicon genome size passed to Canu. Use Canu notation (e.g. `1.5k`, `1500`). |
| `--polishing_reads` | `100` | Number of corrected reads used to build the consensus. |
| `--stopOnLowCoverage` | `1` | Canu `stopOnLowCoverage` setting. |
| `--minInputCoverage` | `2` | Canu `minInputCoverage` setting. |
| `--minReadLength` | `500` | Canu `minReadLength` setting (bp). |
| `--minOverlapLength` | `200` | Canu `minOverlapLength` setting (bp). |
| `--useGrid` | `false` | Canu `useGrid` setting. Set to `false` when running in containerised environments. |

### Classification

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--classification` | `kraken2` | Classification tool(s): `kraken2`, `blast`, `seqmatch`, or `full` (all three in parallel). |
| `--kraken2_db` | `false` | Path to a Kraken2 database directory. Required for `kraken2` and `full` classification, and for `--remove_unclassified`. |
| `--blast_db` | `false` | Path to BLAST database prefix (e.g. `/path/to/db/16S_ribosomal_RNA`). The directory is derived automatically. Required for `blast` and `full`. |
| `--seqmatch_db` | `false` | Path to an RDP SeqMatch-compatible database file. Required for `seqmatch` and `full`. |
| `--seqmatch_accession` | `false` | Path to an RDP accession-to-taxonomy TSV mapping file. Required for `seqmatch` and `full`. |
| `--taxonomy` | `false` | Path to a tab-delimited taxonomy file used by `JOIN_RESULTS` to annotate hits with lineage (name \| species \| genus \| family \| order). |
| `--remove_unclassified` | `false` | Run Kraken2 on QC-passed reads before clustering to remove non-target sequences. Requires `--kraken2_db`. |
| `--reclassifyOnFail` | `false` | If Kraken2 does not classify to species level, automatically retry with SeqMatch. Requires `--seqmatch_db` and `--seqmatch_accession`. |

#### Obtaining classification databases

- **Kraken2**: Download a pre-built database from the [Kraken2 library](https://benlangmead.github.io/aws-indexes/k2) (e.g. `16S_Silva` or `PlusPF`), or build one with `kraken2-build`.
- **BLAST**: Download `16S_ribosomal_RNA` from [NCBI BLAST FTP](https://ftp.ncbi.nlm.nih.gov/blast/db/) (`update_blastdb.pl 16S_ribosomal_RNA`). Supply the full path prefix (without file extension).
- **SeqMatch / RDP**: Download the RDP training set and accession files from the [RDP repository](https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/). The `--seqmatch_accession` file should be tab-separated with columns `accession`, `taxid`, `name`, `species`, `genus`, `family`, `order`.
- **Taxonomy**: The taxonomy file is a tab-delimited file where the first column is `taxid` and subsequent pipe-delimited fields are the lineage. It can be derived from NCBI `nodes.dmp` / `names.dmp` using standard nf-core/taxprofiler or custom parsing scripts.

---

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/nanopath --input samplesheet.csv --outdir <OUTDIR> -profile docker
```

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> ⚠️ Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
> The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/nanopath -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
classification: 'kraken2'
kraken2_db: '/path/to/kraken2_db'
taxonomy: '/path/to/taxonomy.tsv'
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/nanopath
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/nanopath releases page](https://github.com/nf-core/nanopath/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> 💡 If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

> **Note:** This pipeline has been developed and tested with **Docker**, **Singularity**, and **Conda/Mamba** only. Other profiles listed below are included as standard nf-core defaults but have **not been tested** with nf-core/nanopath and may not work as expected.

Note that multiple profiles can be loaded, for example: `-profile test,docker` — the order of arguments is important; they are loaded in sequence so later profiles can overwrite earlier ones.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

**Tested and supported:**

- `docker`
  - Run with [Docker](https://docker.com/). Recommended for most users.
- `singularity`
  - Run with [Singularity](https://sylabs.io/docs/). Recommended for HPC environments.
- `conda` / `mamba`
  - Run with [Conda](https://conda.io/docs/) or [Mamba](https://mamba.readthedocs.io/). Use only when Docker or Singularity are not available.
- `test`
  - A minimal test profile for automated testing.

**Not tested with this pipeline (standard nf-core defaults):**

- `podman` — [Podman](https://podman.io/)
- `shifter` — [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud` — [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer` — [Apptainer](https://apptainer.org/)

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
