#!/bin/bash
#SBATCH --job-name=rygb_dl_cocktail_resume
#SBATCH --partition=Orion
#SBATCH --array=1-6
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/%u/rygb/logs/download_cocktail_resume_%A_%a.out
#SBATCH --error=/scratch/%u/rygb/logs/download_cocktail_resume_%A_%a.out
#
# One-off resume/parallelization for COCKTAIL_GBP's download. 02b_download_wgs.sbatch's
# task 5 (408 runs, one serial task) hit its 24h time limit at 112/408 runs (job
# 27788909, cancelled 2026-11-15T18:48 DUE TO TIME LIMIT). At the observed rate
# (~4.7 runs/hr), finishing the remaining ~296 serially would take another 2-3 days.
#
# This splits the FULL run list into 6 parallel chunks (by position modulo 6) and
# reuses 02b's exact per-run download logic (prefetch + fasterq-dump + pigz), so it's
# safely resumable: any run whose _R1/_R2.fastq.gz already exists (from 02b's earlier
# partial run) is skipped, not re-downloaded. Submit once; delete once COCKTAIL_GBP's
# raw/ pair count reaches 408 (see 02b_download_wgs.sbatch's own EXPECT=408).

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
module load sra-tools/3.1.0

STUDY=COCKTAIL_GBP
N_TASKS=6
OUT=$DATA_ROOT/raw/$STUDY
LIST=$OUT/runs.tsv   # same master list 02b_download_wgs.sbatch already built
TMP=$DATA_ROOT/work/dl_${STUDY}_resume_${SLURM_ARRAY_TASK_ID}
mkdir -p "$OUT" "$TMP"
cd "$TMP" || exit 1
THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "=== $STUDY resume chunk $SLURM_ARRAY_TASK_ID/$N_TASKS   $(date)"
have0=$(ls "$OUT"/*_R1.fastq.gz 2>/dev/null | wc -l)
echo "already present (any chunk, from the earlier partial run): $have0 / 408"

fail=0; done_n=0; skip=0; idx=0
while IFS=$'\t' read -r run; do
  [ -z "$run" ] && continue
  idx=$((idx+1))
  [ $(( (idx-1) % N_TASKS + 1 )) -eq "$SLURM_ARRAY_TASK_ID" ] || continue

  target1=$OUT/${run}_R1.fastq.gz
  target2=$OUT/${run}_R2.fastq.gz
  # gzip -t, not just non-empty — see 02b_download_wgs.sbatch's comment on this same check.
  if [ -s "$target1" ] && [ -s "$target2" ] && gzip -t "$target1" 2>/dev/null && gzip -t "$target2" 2>/dev/null; then
    skip=$((skip+1)); continue
  fi

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
    echo "FAIL not paired-end (or missing reads) for $run"
    fail=$((fail+1))
  fi
  rm -rf "$TMP/$run" "$TMP/${run}"_*.fastq
done < "$LIST"

echo
echo "=== chunk $SLURM_ARRAY_TASK_ID done  $(date)"
echo "downloaded now: $done_n | already present: $skip | failures: $fail"
have=$(ls "$OUT"/*_R1.fastq.gz 2>/dev/null | wc -l)
echo "total pairs in $OUT now: $have (expect 408)"
rm -rf "$TMP"
[ "$fail" -eq 0 ] || exit 1
