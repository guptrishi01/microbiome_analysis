#!/bin/bash
#SBATCH --job-name=rygb_dada2
#SBATCH --partition=Orion
#SBATCH --array=1-4
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/%u/rygb/logs/dada2_%A_%a.out
#SBATCH --error=/scratch/%u/rygb/logs/dada2_%A_%a.out
#
# DADA2 on forward reads, one array task per 16S study.
# Per-study trimming derived from the repo's own ASV tables (uniform lengths 200/200/180/150):
#   BS16S  cutadapt through 515F, truncLen 200            -> 200 bp ASVs
#   Assal  cutadapt through 515F, truncLen 200            -> 200 bp ASVs
#   Afshar trimLeft 20, truncLen 200                      -> 180 bp ASVs
#   Ilhan  no trim, truncLen 150                          -> 150 bp ASVs
# Produces ForwardReads.txt + taxForwardReads.txt in the repo's layout, keeps the
# filtered FASTQs (Kraken2 reads those next), and reports overlap with the repo's ASVs.

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
ENV16S=$SW_ROOT/envs/rygb-16s3
REPO=$DATA_ROOT/RYGB_IntegratedAnalysis2020
THREADS=${SLURM_CPUS_PER_TASK:-16}

case $SLURM_ARRAY_TASK_ID in
  1) STUDY=BS16S;  REPODIR=RYGB_BS;         PRIMER=GTGYCAGCMGCCGCGGTAA; TRUNC=200; TRIMLEFT=0;  EXPECT_LEN=200; REPO_ASVS=1525 ;;
  2) STUDY=Assal;  REPODIR=RYGB_Assal2020;  PRIMER=GTGYCAGCMGCCGCGGTAA; TRUNC=200; TRIMLEFT=0;  EXPECT_LEN=200; REPO_ASVS=3274 ;;
  3) STUDY=Afshar; REPODIR=RYGB_Afshar2018; PRIMER=NONE;                TRUNC=200; TRIMLEFT=20; EXPECT_LEN=180; REPO_ASVS=3613 ;;
  4) STUDY=Ilhan;  REPODIR=RYGB_Ilhan2020;  PRIMER=NONE;                TRUNC=150; TRIMLEFT=0;  EXPECT_LEN=150; REPO_ASVS=5378 ;;
  *) echo "bad array index"; exit 1 ;;
esac

RAW=$DATA_ROOT/raw/$STUDY
CUT=$DATA_ROOT/work/cutadapt/$STUDY
FILT=$DATA_ROOT/work/filtered/$STUDY
OUT=$DATA_ROOT/input/$REPODIR/DADA2
mkdir -p "$CUT" "$FILT" "$OUT"
echo "=== $STUDY  primer=$PRIMER truncLen=$TRUNC trimLeft=$TRIMLEFT  expect ${EXPECT_LEN}bp ASVs   $(date)"
echo "raw files: $(ls "$RAW"/*_R1.fastq.gz | wc -l)"

# ---------------------------------------------------------- primer position check
echo
echo "--- where does the 16S conserved region sit in the raw reads? (first file, 5000 reads)"
zcat "$(ls "$RAW"/*_R1.fastq.gz | head -1)" | awk 'NR%4==2' | head -5000 \
  | awk '{i=index($0,"GCCGCGGTAA"); print (i? i : "no-515F-anchor")}' \
  | sort | uniq -c | sort -rn | head -4 | sed 's/^/    /'

# ---------------------------------------------------------- cutadapt (BS, Assal)
if [ "$PRIMER" != "NONE" ]; then
  echo
  echo "--- cutadapt: remove everything up to and including $PRIMER"
  n=0
  for f in "$RAW"/*_R1.fastq.gz; do
    s=$(basename "$f" _R1.fastq.gz)
    [ -s "$CUT/${s}_R1.fastq.gz" ] && { n=$((n+1)); continue; }
    "$ENV16S/bin/cutadapt" -g "$PRIMER" --discard-untrimmed -j "$THREADS" \
      -o "$CUT/${s}_R1.fastq.gz" "$f" > "$CUT/${s}.cutadapt.log" 2>&1 || { echo "FAIL cutadapt $s"; exit 1; }
    n=$((n+1))
  done
  echo "    trimmed $n samples"
  echo "    retention across samples (reads written / reads processed):"
  grep -h "Reads with adapters\|Reads written" "$CUT"/*.cutadapt.log | head -4 | sed 's/^/      /'
  grep -h "Reads written" "$CUT"/*.cutadapt.log | sed 's/.*(\([0-9.]*\)%).*/\1/' \
    | sort -n | awk '{a[NR]=$1} END{print "      min "a[1]"%  median "a[int(NR/2)]"%  max "a[NR]"%"}'
  INDIR=$CUT
else
  INDIR=$RAW
fi

# ---------------------------------------------------------- DADA2
cat > "$DATA_ROOT/work/dada2_$STUDY.R" <<'RSCRIPT'
suppressMessages(library(dada2)); suppressMessages(library(ShortRead))
a <- commandArgs(trailingOnly=TRUE)
inDir <- a[1]; filtDir <- a[2]; outDir <- a[3]
truncLen <- as.integer(a[4]); trimLeft <- as.integer(a[5])
silvaTrain <- a[6]; silvaSpecies <- a[7]; threads <- as.integer(a[8])
repoTable <- a[9]; expectLen <- as.integer(a[10])

fns <- sort(list.files(inDir, pattern="_R1.fastq.gz$", full.names=TRUE))
samples <- sub("_R1.fastq.gz$", "", basename(fns))
filts <- file.path(filtDir, paste0(samples, "_R1.fastq.gz"))
cat("samples:", length(samples), "\n")

# paper: filterAndTrim with defaults and maxEE = 2
track <- filterAndTrim(fns, filts, truncLen=truncLen, trimLeft=trimLeft,
                       maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                       compress=TRUE, multithread=FALSE, verbose=FALSE)
rownames(track) <- samples
cat("reads in/out — median retention:",
    round(median(track[,2]/track[,1])*100,1), "%\n")
lost <- rownames(track)[track[,2] < 1000]
if (length(lost)) cat("samples under 1000 reads after filtering (norm() drops these):",
                      paste(lost, collapse=", "), "\n")

keep <- file.exists(filts)
filts <- filts[keep]; samples <- samples[keep]

err <- learnErrors(filts, multithread=threads, verbose=FALSE)
dd  <- dada(filts, err=err, multithread=threads, verbose=FALSE)
seqtab <- makeSequenceTable(dd)
rownames(seqtab) <- samples
cat("ASVs before chimera removal:", ncol(seqtab), "\n")

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=threads, verbose=FALSE)
cat("ASVs after chimera removal:", ncol(seqtab.nochim),
    "| reads kept:", round(sum(seqtab.nochim)/sum(seqtab)*100,1), "%\n")

lens <- table(nchar(colnames(seqtab.nochim)))
cat("ASV length distribution:\n"); print(lens)
if (!all(nchar(colnames(seqtab.nochim)) == expectLen))
  cat("!! WARNING: expected every ASV to be", expectLen, "bp\n")

# taxonomy — SILVA 132, as in the paper
tax <- assignTaxonomy(colnames(seqtab.nochim), silvaTrain, multithread=threads, verbose=FALSE)
tax <- addSpecies(tax, silvaSpecies, verbose=FALSE)
rownames(tax) <- colnames(seqtab.nochim)

# repo layout: write.table defaults reproduce the committed files exactly
write.table(seqtab.nochim, file.path(outDir, "ForwardReads.txt"), sep="\t", quote=TRUE)
write.table(tax,           file.path(outDir, "taxForwardReads.txt"), sep="\t", quote=TRUE)
write.table(track,         file.path(outDir, "filter_track.txt"), sep="\t", quote=FALSE)

# how close are we to the published ASV set?
if (file.exists(repoTable)) {
  old <- colnames(read.table(repoTable, sep="\t", header=TRUE, nrows=1, check.names=FALSE))
  new <- colnames(seqtab.nochim)
  shared <- intersect(toupper(old), toupper(new))
  cat("\nrepo ASVs:", length(old), "| ours:", length(new),
      "| shared:", length(shared),
      sprintf("(%.1f%% of repo, %.1f%% of ours)\n",
              100*length(shared)/length(old), 100*length(shared)/length(new)))
}
cat("\nwrote:", list.files(outDir), "\n")
RSCRIPT

echo
echo "--- DADA2"
"$ENV16S/bin/Rscript" "$DATA_ROOT/work/dada2_$STUDY.R" \
  "$INDIR" "$FILT" "$OUT" "$TRUNC" "$TRIMLEFT" \
  "$SW_ROOT/db/silva132/silva_nr_v132_train_set.fa.gz" \
  "$SW_ROOT/db/silva132/silva_species_assignment_v132.fa.gz" \
  "$THREADS" "$REPO/input/$REPODIR/DADA2/ForwardReads.txt" "$EXPECT_LEN"
rc=$?

echo
echo "=== $STUDY done (rc=$rc)  $(date)"
echo "repo had $REPO_ASVS ASVs for this study"
ls -lh "$OUT" | sed 's/^/   /'
echo "filtered reads kept for Kraken2: $(ls "$FILT"/*.fastq.gz 2>/dev/null | wc -l) files, $(du -sh "$FILT" | cut -f1)"
exit $rc
