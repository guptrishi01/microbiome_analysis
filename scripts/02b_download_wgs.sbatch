#!/bin/bash
#SBATCH --job-name=rygb_dl_wgs
#SBATCH --partition=Orion
#SBATCH --array=1-6
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/%u/rygb/logs/download_wgs_%A_%a.out
#SBATCH --error=/scratch/%u/rygb/logs/download_wgs_%A_%a.out
#
# WGS taxonomy-track downloads: the 4 NEW shotgun-metagenomic cohorts identified this
# session (see MACRIB_Research_Statement.docx's "Planned Future Work" section and
# cohorts/<Study>/NEEDS_METADATA.md). Unlike 02_download.sbatch (16S, forward-reads-only,
# matching the original authors' DADA2-on-forward-reads convention), this keeps BOTH
# paired-end reads — shotgun classification (Kraken2/Bracken) and kneaddata both expect
# paired input, and discarding R2 would throw away half the sequencing depth for no reason.
#
#   task 1 Penney_RYGB  PRJNA473348    25 runs   Pre/POST/1YrPOST; Arm=RYGB — PARTIALLY
#                                                 RESOLVED 2026-11-14 via the paper's own
#                                                 Mendeley deposit (10.17632/t76nm3yfzh.3);
#                                                 only 10 of the paper's 13 reported RYGB
#                                                 patients recovered — see cohorts/
#                                                 Penney_RYGB/NEEDS_METADATA.md. The other
#                                                 ~64 of the original 120 runs (unresolved
#                                                 patients) are intentionally NOT downloaded.
#   task 2 Penney_SG    PRJNA473348    31 runs   Pre/POST/1YrPOST; Arm=SG — RESOLVED, matches
#                                                 the paper's reported SG arm exactly (14/14).
#                                                 Same source/caveat as Penney_RYGB above.
#   task 3 Oseberg_RYGB PRJEB59861     78 runs   B/12M; Arm=RYGB — RESOLVED 2026-11-14 via
#                                                 the paper's own Zenodo source-data deposit
#                                                 (10.5281/zenodo.19203868); cohorts/Oseberg/
#                                                 was split into Oseberg_RYGB/Oseberg_SG.
#   task 4 Oseberg_SG   PRJEB59861     76 runs   B/12M; Arm=SG — RESOLVED, same source as above.
#   task 5 COCKTAIL_GBP PRJNA904651   408 runs   time/ID/Arm ALL unknown (anonymized
#                                                 aliases) — see cohorts/COCKTAIL_GBP/
#                                                 NEEDS_METADATA.md. NOTE: 408 runs is the
#                                                 FULL BioProject (GBP + excluded VLED diet
#                                                 arm, since the two arms aren't
#                                                 distinguishable from public metadata
#                                                 either) — do not assume every run here is
#                                                 a surgical sample until that's resolved.
#   task 6 Hayashi_LSG  PRJNA1329344   38 runs   time/ID resolved (BioSample metadata);
#                                                 Arm (LSG/IGBP) still unknown — see
#                                                 cohorts/Hayashi_LSG/NEEDS_METADATA.md.
#                                                 NOTE: 38 runs is the FULL BioProject
#                                                 (LSG + excluded IGBP arm, same caveat as
#                                                 COCKTAIL above).
#
# BS-WGS (PRJNA668357/PRJNA668472, the Fouladi et al. 2021 metagenomic subset of the
# existing phase-1 BS cohort) is deliberately NOT downloaded here: the original authors'
# repo already ships fully-processed Kraken2 species/genus/etc. count tables and curated
# metadata for it at $DATA_ROOT/RYGB_IntegratedAnalysis2020/input/RYGB_BS_Metagenomics/ —
# same "reuse the authors' own tables" treatment phase-1's 16S cohorts already get. See
# analysis/RScripts/compareStudies/compareStudies_Metagenomics.R (this repo's generalized
# version) for how it's wired into the cross-cohort comparison without reprocessing.
#
# Resumable: a run whose _R1/_R2 fastq.gz already exist is skipped.
# Submit:  sbatch 02b_download_wgs.sbatch

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
COHORTS=$SW_ROOT/cohorts
module load sra-tools/3.1.0

case $SLURM_ARRAY_TASK_ID in
  1) STUDY=Penney_RYGB;  EXPECT=25  ;;
  2) STUDY=Penney_SG;    EXPECT=31  ;;
  3) STUDY=Oseberg_RYGB; EXPECT=78  ;;
  4) STUDY=Oseberg_SG;   EXPECT=76  ;;
  5) STUDY=COCKTAIL_GBP; EXPECT=408 ;;
  6) STUDY=Hayashi_LSG;  EXPECT=38  ;;
  *) echo "bad array index"; exit 1 ;;
esac

OUT=$DATA_ROOT/raw/$STUDY
TMP=$DATA_ROOT/work/dl_$STUDY
THREADS=${SLURM_CPUS_PER_TASK:-8}
mkdir -p "$OUT" "$TMP"
cd "$TMP" || exit 1
echo "=== $STUDY  expecting $EXPECT runs   $(date)"

# run list: column 1 = Run, from this repo's own cohorts/<Study>/metaData.txt
LIST=$OUT/runs.tsv
if [ -s "$LIST" ]; then
  echo "run list exists: $LIST"
else
  tail -n +2 "$COHORTS/$STUDY/metaData.txt" | cut -f1 | grep . > "$LIST"
fi

N=$(wc -l < "$LIST")
echo "run list: $N entries (expected $EXPECT)"
[ "$N" -eq "$EXPECT" ] || echo "!! WARNING: count mismatch — check $LIST before trusting downstream output"
head -3 "$LIST"

# ---------------------------------------------------------------- download (paired-end)
fail=0; done_n=0; skip=0
while IFS=$'\t' read -r run; do
  [ -z "$run" ] && continue
  target1=$OUT/${run}_R1.fastq.gz
  target2=$OUT/${run}_R2.fastq.gz
  if [ -s "$target1" ] && [ -s "$target2" ]; then skip=$((skip+1)); continue; fi

  if ! prefetch --max-size 100G -O "$TMP" "$run" >"$TMP/$run.prefetch.log" 2>&1; then
    echo "FAIL prefetch $run"; tail -3 "$TMP/$run.prefetch.log"; fail=$((fail+1)); continue
  fi
  if ! fasterq-dump --split-files --threads 8 -O "$TMP" "$TMP/$run" >"$TMP/$run.fqdump.log" 2>&1; then
    echo "FAIL fasterq-dump $run"; tail -3 "$TMP/$run.fqdump.log"; fail=$((fail+1)); rm -rf "$TMP/$run"; continue
  fi

  fwd=$TMP/${run}_1.fastq
  rev=$TMP/${run}_2.fastq
  if [ -s "$fwd" ] && [ -s "$rev" ]; then
    "$SW_ROOT/envs/rygb-16s/bin/pigz" -p "$THREADS" -c "$fwd" > "$target1"
    "$SW_ROOT/envs/rygb-16s/bin/pigz" -p "$THREADS" -c "$rev" > "$target2"
    done_n=$((done_n+1))
  else
    echo "FAIL not paired-end (or missing reads) for $run — single-end submissions aren't"
    echo "     expected for any of these 4 cohorts; investigate before assuming it's fine"
    fail=$((fail+1))
  fi
  rm -rf "$TMP/$run" "$TMP/${run}"_*.fastq
done < "$LIST"

# ---------------------------------------------------------------- report
echo
echo "=== $STUDY done  $(date)"
echo "downloaded now: $done_n | already present: $skip | failures: $fail"
have=$(ls "$OUT"/*_R1.fastq.gz 2>/dev/null | wc -l)
echo "pairs in $OUT: $have (expected $EXPECT)"
du -sh "$OUT"
if [ "$have" -gt 0 ]; then
  f=$(ls "$OUT"/*_R1.fastq.gz | head -1)
  echo "read length check ($(basename "$f")):"
  zcat "$f" | head -2 | tail -1 | awk '{print "   first read: "length($0)" bp"}'
fi
rm -rf "$TMP"
[ "$fail" -eq 0 ] || exit 1
