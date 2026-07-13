#!/bin/bash
#SBATCH --job-name=db_kraken2_16s
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=db_kraken2_16s.%j.out
#SBATCH --error=db_kraken2_16s.%j.err
#SBATCH --mail-user=j.goodall@sheffield.ac.uk
#SBATCH --mail-type=END,FAIL

# Build the kraken2 16S index from the SAME curated NCBI 16S_ribosomal_RNA set that
# BLAST uses (bacteria+archaea 16S; no 18S/fungi), with NCBI taxonomy so taxids
# resolve in rankedlineage.dmp. Tools (kraken2 + blastdbcmd) come from the pipeline's
# apptainer image -- no conda.
set -euo pipefail

DEST=/mnt/parscratch/users/md1jgzz/fresh_databases
DATE=$(date +%Y%m%d)
DB="${DEST}/kraken2_16S_${DATE}"
THREADS=${SLURM_CPUS_PER_TASK:-8}
IMG=/mnt/parscratch/users/md1jgzz/singularity_cache/mbdabrowska1-full-classification-1.0.img

mkdir -p "$DEST" "$(dirname "$IMG")"

# --- container runtime (apptainer or singularity) ---
command -v apptainer >/dev/null 2>&1 || module load apptainer 2>/dev/null || true
RUNTIME=$(command -v apptainer || command -v singularity || true)
[ -n "$RUNTIME" ] || { echo "ERROR: no apptainer/singularity on PATH"; exit 1; }
# pull the image if it isn't cached yet
[ -f "$IMG" ] || "$RUNTIME" pull "$IMG" docker://mbdabrowska1/full-classification:1.0
# run a command inside the image, with parscratch bound in
run() { "$RUNTIME" exec -B /mnt/parscratch "$IMG" "$@"; }

# --- 1) NCBI 16S_ribosomal_RNA (host wget/tar) ---
WORK="${DB}/build"; mkdir -p "$WORK"; cd "$WORK"
wget -q https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
tar xzf 16S_ribosomal_RNA.tar.gz && rm -f 16S_ribosomal_RNA.tar.gz

# --- 2) dump sequences with explicit NCBI taxids in kraken2 header format ---
run blastdbcmd -db 16S_ribosomal_RNA -entry all -outfmt '%a %T %s' \
  | awk '{print ">"$1"|kraken:taxid|"$2"\n"$3}' > 16S.fasta

# --- 3) build with NCBI taxonomy (--skip-maps: taxids are already in the headers) ---
cd "$DB"
run kraken2-build --db "$DB" --download-taxonomy --skip-maps --threads "$THREADS"
run kraken2-build --db "$DB" --add-to-library "${WORK}/16S.fasta"
run kraken2-build --db "$DB" --build --threads "$THREADS"
run kraken2-build --db "$DB" --clean
rm -rf "$WORK"

echo "Done."
echo "  --kraken2_db ${DB}"
echo "Sanity:  ${RUNTIME} exec ${IMG} kraken2-inspect --db ${DB} | head"
