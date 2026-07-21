#!/bin/bash
# Run nanopath on the HPC with compute capped to match the GridION (GXB03144), so an HPC
# result can be compared fairly against the on-instrument result.
#
# GridION specs read off the device:
#   CPU : Intel Core i7-10700K -- 8 cores / 16 threads  (nproc = 16)
#   RAM : 62 GiB total
#   GPU : Quadro GV100 (32 GB)  -- BASECALLING ONLY, not used by nanopath's pipeline
#
# So we request 16 CPUs + 62G, and cap Nextflow to the same via --max_cpus/--max_memory
# (a big SLURM allocation alone won't stop Nextflow scaling past the GridION).

#SBATCH --job-name=nanopath_gridion_equiv
#SBATCH --cpus-per-task=16
#SBATCH --mem=62G
#SBATCH --time=12:00:00
#SBATCH --output=nanopath_gridion_equiv.%j.out
#SBATCH --error=nanopath_gridion_equiv.%j.err
#SBATCH --mail-user=j.goodall@sheffield.ac.uk
#SBATCH --mail-type=END,FAIL
## No GPU: nanopath's clustering/consensus/classification are CPU-only. If you ever ALSO
## basecall on the HPC, add:  #SBATCH --gres=gpu:1   (note: HPC GPUs != the GridION's GV100)

set -euo pipefail

# ---- environment (match how you normally launch nanopath) ----
module load Nextflow/25.10.2
export NXF_SINGULARITY_CACHEDIR=/mnt/parscratch/users/md1jgzz/singularity_cache
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# ---- EDIT THESE ----------------------------------------------------------------
PIPELINE=/path/to/nf-core-nanopath                 # repo dir (or "." if you sbatch from it)
INPUT=/path/to/experiment_info.xlsx                # samplesheet
OUTDIR=/mnt/parscratch/users/md1jgzz/gridion_equiv_out

# databases -- default to the NEW dense genome-derived dbs; point at the old sparse ones to A/B
KRAKEN2_DB=/mnt/parscratch/users/md1jgzz/fresh_databases/full/kraken2_genome16S_20260714
BLAST_DB=/mnt/parscratch/users/md1jgzz/fresh_databases/full/genome_16S_blast_20260714/genome_16S
SEQMATCH_DB=/mnt/parscratch/users/md1jgzz/fresh_databases/CHANGE_ME_seqmatch_db
SEQMATCH_ACC=/mnt/parscratch/users/md1jgzz/fresh_databases/CHANGE_ME_accession.txt
TAXONOMY=/mnt/parscratch/users/md1jgzz/fresh_databases/CHANGE_ME_rankedlineage.dmp
# -------------------------------------------------------------------------------

cd "$PIPELINE"

# --max_cpus/--max_memory are the GridION cap. --mem=62G but cap processes at 60G so the
# Nextflow driver + overhead still fit inside the allocation.
nextflow run . -profile singularity -resume \
  --max_cpus 16 --max_memory 60.GB \
  --classification full \
  --input        "$INPUT" \
  --kraken2_db   "$KRAKEN2_DB" \
  --blast_db     "$BLAST_DB" \
  --seqmatch_db  "$SEQMATCH_DB" \
  --seqmatch_accession "$SEQMATCH_ACC" \
  --taxonomy     "$TAXONOMY" \
  --outdir       "$OUTDIR"

echo "Done. Results in $OUTDIR"
