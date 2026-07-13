#!/bin/bash
#SBATCH --job-name=db_kraken2_16s
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=db_kraken2_16s.%j.out
#SBATCH --error=db_kraken2_16s.%j.err
#SBATCH --mail-user=j.goodall@sheffield.ac.uk
#SBATCH --mail-type=END,FAIL

set -euo pipefail

DEST=/mnt/parscratch/users/md1jgzz/fresh_databases
DATE=$(date +%Y%m%d)
DB="${DEST}/kraken2_16S_${DATE}"
THREADS=${SLURM_CPUS_PER_TASK:-8}

# tools: kraken2 (+ blast for blastdbcmd/dustmasker)
eval "$(conda shell.bash hook)"
[ -d "${DEST}/.dbtools" ] || conda create -y -p "${DEST}/.dbtools" -c conda-forge -c bioconda kraken2 blast wget
source activate "${DEST}/.dbtools"

WORK="${DB}/build"; mkdir -p "$WORK"; cd "$WORK"

# 1) get NCBI 16S_ribosomal_RNA (same product BLAST uses)
wget -q https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
tar xzf 16S_ribosomal_RNA.tar.gz && rm -f 16S_ribosomal_RNA.tar.gz

# 2) dump the sequences with explicit NCBI taxids in kraken2 header format
blastdbcmd -db 16S_ribosomal_RNA -entry all -outfmt '%a %T %s' \
  | awk '{print ">"$1"|kraken:taxid|"$2"\n"$3}' > 16S.fasta

# 3) build with NCBI taxonomy; --skip-maps because taxids are already in the headers
cd "$DB"
kraken2-build --db "$DB" --download-taxonomy --skip-maps --threads "$THREADS"
kraken2-build --db "$DB" --add-to-library "${WORK}/16S.fasta"
kraken2-build --db "$DB" --build --threads "$THREADS"
kraken2-build --db "$DB" --clean
rm -rf "$WORK"

echo "Done."
echo "  --kraken2_db ${DB}"
echo "Sanity check:  kraken2-inspect --db ${DB} | head"
