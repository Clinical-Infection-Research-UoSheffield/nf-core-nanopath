#!/bin/bash
#SBATCH --job-name=db_genome16s
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=db_genome16s.%j.out
#SBATCH --error=db_genome16s.%j.err
#SBATCH --mail-user=j.goodall@sheffield.ac.uk
#SBATCH --mail-type=END,FAIL

# Build a DENSE 16S reference from RefSeq genomes (not the sparse Targeted Loci db).
# Rationale: Targeted Loci carries ~1 curated 16S per species; a real clinical strain
# whose 16S drifts from that lone record slides onto a neighbour species (the Hflu case).
# RefSeq COMPLETE genomes give thousands of strains for common pathogens x ~6 operons each,
# so exact-identity matches actually exist. We dereplicate to keep every UNIQUE 16S variant
# and drop only exact duplicate copies -> dense where it matters, still slim.
#
# The single output FASTA feeds BOTH makeblastdb and kraken2-build, so BLAST and kraken2
# stop sharing the type-strain blind spot. Taxids are native NCBI (no SILVA remap).

set -euo pipefail

DEST=/mnt/parscratch/users/md1jgzz/fresh_databases
OUT="${DEST}/full"          # keep these dense "full-genome-derived" dbs separate from the others
mkdir -p "$OUT"
DATE=$(date +%Y%m%d)

# assembly quality to include: "Complete Genome" alone is leaner; add "Chromosome" for more depth
LEVELS='Complete Genome|Chromosome'
MIN_LEN=1200            # drop partial 16S fragments (full-length ~1500 bp)
DOWNLOAD_JOBS=12        # parallel wget streams -- keep modest to stay polite to NCBI
THREADS=${SLURM_CPUS_PER_TASK:-16}

WORK="${OUT}/genome16s_build_${DATE}"
PARTS="${WORK}/parts"
mkdir -p "$PARTS"
cd "$WORK"

# ---- tools: use the conda env's binaries directly (no `conda` command in the batch job) ----
DBTOOLS="${DEST}/.dbtools"
[ -x "${DBTOOLS}/bin/kraken2-build" ] || { echo "ERROR: conda env missing at ${DBTOOLS}; create it interactively first: conda create -y -p ${DBTOOLS} -c conda-forge -c bioconda kraken2 blast rsync"; exit 1; }
export PATH="${DBTOOLS}/bin:${PATH}"

# ---- 1) list target assemblies from the RefSeq bacteria summary ----
echo "[1/5] fetching assembly summary"
wget -q -O assembly_summary.txt https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

# cols: 6=taxid, 11=version_status, 12=assembly_level, 20=ftp_path. Filter on assembly_level
# (a QUALITY filter that keeps ALL strains) -- NOT refseq_category, which is 1-per-species.
awk -F'\t' -v levels="$LEVELS" '
  /^#/ { next }
  $11=="latest" && $20!="na" && ($12 ~ "^(" levels ")$") { print $6, $20 }
' assembly_summary.txt > targets.txt
echo "      $(wc -l < targets.txt) assemblies to scan"

# ---- 2) pull 16S from each assembly's annotated rna_from_genomic file, tag with its taxid ----
# We fetch only the small *_rna_from_genomic.fna.gz (annotated rRNA/tRNA/ncRNA), not whole genomes,
# and keep records whose defline says [product=16S ribosomal RNA]. Output: <taxid>\t<linearised seq>.
echo "[2/5] downloading + extracting 16S (${DOWNLOAD_JOBS} parallel)"
export PARTS MIN_LEN

fetch_one() {
  local taxid="$1" ftp="$2" base url tmp
  base=$(basename "$ftp")
  url="${ftp}/${base}_rna_from_genomic.fna.gz"
  tmp=$(mktemp)
  if wget -q --tries=2 -O "${tmp}.gz" "$url" 2>/dev/null && [ -s "${tmp}.gz" ]; then
    gzip -dc "${tmp}.gz" 2>/dev/null | awk -v taxid="$taxid" -v minlen="$MIN_LEN" '
      /^>/ {
        if (keep && length(seq) >= minlen) print taxid "\t" seq
        seq=""; keep = ($0 ~ /\[product=16S ribosomal RNA\]/); next
      }
      { seq = seq $0 }
      END { if (keep && length(seq) >= minlen) print taxid "\t" seq }
    ' >> "${PARTS}/part.${BASHPID}.tsv"
  fi
  rm -f "$tmp" "${tmp}.gz"
  return 0
}
export -f fetch_one

xargs -a targets.txt -P "$DOWNLOAD_JOBS" -n 2 bash -c 'fetch_one "$1" "$2"' _ || true

cat "${PARTS}"/part.*.tsv > raw.tsv 2>/dev/null || { echo "ERROR: no 16S extracted"; exit 1; }
echo "      $(wc -l < raw.tsv) raw 16S sequences"

# ---- 3) dereplicate: unique (taxid, sequence) rows ----
# Collapses the many identical operon copies / identical strains, but KEEPS a sequence that
# occurs under different taxids as separate rows (kraken2 resolves shared k-mers to the LCA
# at build time; BLAST just reports each). This is the "slim without going sparse" step.
echo "[3/5] dereplicating"
LC_ALL=C sort -u --parallel="$THREADS" -S 50% raw.tsv > derep.tsv
echo "      $(wc -l < derep.tsv) unique 16S variants"

# ---- 4) emit the build inputs from the derep'd set ----
# blast.fasta + taxid_map.txt for makeblastdb; kraken.fasta (headers carry the taxid) for kraken2.
rm -f blast.fasta taxid_map.txt kraken.fasta
awk -F'\t' '{
  id="S" NR; tax=$1; seq=$2
  print ">" id "\n" seq                     >> "blast.fasta"
  print id "\t" tax                          >> "taxid_map.txt"
  print ">" id "|kraken:taxid|" tax "\n" seq >> "kraken.fasta"
}' derep.tsv

# ---- 5a) BLAST db ----
echo "[4/5] building BLAST db"
BLAST_OUT="${OUT}/genome_16S_blast_${DATE}"
mkdir -p "$BLAST_OUT"
makeblastdb -in blast.fasta -dbtype nucl -parse_seqids \
  -taxid_map taxid_map.txt -title "genome_16S_${DATE}" \
  -out "${BLAST_OUT}/genome_16S"

# ---- 5b) kraken2 db (NCBI taxonomy; taxids already in the headers -> --skip-maps) ----
echo "[5/5] building kraken2 db"
KDB="${OUT}/kraken2_genome16S_${DATE}"
mkdir -p "$KDB"
kraken2-build --db "$KDB" --download-taxonomy --skip-maps --use-ftp --threads "$THREADS"
kraken2-build --db "$KDB" --add-to-library "${WORK}/kraken.fasta"
kraken2-build --db "$KDB" --build --threads "$THREADS"
kraken2-build --db "$KDB" --clean

echo
echo "Done."
echo "  --blast_db   ${BLAST_OUT}/genome_16S"
echo "  --kraken2_db ${KDB}"
echo "Sanity checks:"
echo "  blastdbcmd -db ${BLAST_OUT}/genome_16S -taxids 727 -outfmt '%a %T' | wc -l   # Hflu 16S count"
echo "  kraken2-inspect --db ${KDB} | head"
echo "(remove the build scratch when happy:  rm -rf ${WORK} )"
