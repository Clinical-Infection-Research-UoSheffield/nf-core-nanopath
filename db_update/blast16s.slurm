#!/bin/bash
#SBATCH --job-name=db_blast16s
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=db_blast16s.%j.out
#SBATCH --error=db_blast16s.%j.err
#SBATCH --mail-user=j.goodall@sheffield.ac.uk
#SBATCH --mail-type=END,FAIL

set -euo pipefail

DEST=/mnt/parscratch/users/md1jgzz/fresh_databases
DATE=$(date +%Y%m%d)
DIR="${DEST}/16S_ribosomal_RNA_blast_${DATE}"
mkdir -p "$DIR" && cd "$DIR"

BASE=https://ftp.ncbi.nlm.nih.gov/blast/db
wget -q "${BASE}/16S_ribosomal_RNA.tar.gz"
wget -q "${BASE}/16S_ribosomal_RNA.tar.gz.md5"
md5sum -c 16S_ribosomal_RNA.tar.gz.md5
tar xzf 16S_ribosomal_RNA.tar.gz
rm -f 16S_ribosomal_RNA.tar.gz 16S_ribosomal_RNA.tar.gz.md5

echo "Done."
echo "  --blast_db ${DIR}/16S_ribosomal_RNA"
