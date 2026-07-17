#!/bin/bash
#SBATCH --job-name=db_taxonomy
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=db_taxonomy.%j.out
#SBATCH --error=db_taxonomy.%j.err
#SBATCH --mail-user=j.goodall@sheffield.ac.uk
#SBATCH --mail-type=END,FAIL

set -euo pipefail

DEST=/mnt/parscratch/users/md1jgzz/fresh_databases
DATE=$(date +%Y%m%d)
cd "$DEST"

BASE=https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump
wget -q "${BASE}/new_taxdump.tar.gz"
wget -q "${BASE}/new_taxdump.tar.gz.md5"
md5sum -c new_taxdump.tar.gz.md5
tar xzf new_taxdump.tar.gz rankedlineage.dmp
mv rankedlineage.dmp "rankedlineage_${DATE}.dmp"
rm -f new_taxdump.tar.gz new_taxdump.tar.gz.md5

echo "Done."
echo "  --taxonomy ${DEST}/rankedlineage_${DATE}.dmp"
