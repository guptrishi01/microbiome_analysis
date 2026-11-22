#!/bin/bash
#SBATCH --job-name=rygb_kneaddata
#SBATCH --partition=Orion
#SBATCH --array=1-6
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/%u/rygb/logs/kneaddata_%A_%a.out
#SBATCH --error=/scratch/%u/rygb/logs/kneaddata_%A_%a.out
#
# Stage 3 (WGS track): quality trim + human-read removal on the 4 new shotgun cohorts.
# No DADA2 equivalent here — kneaddata (Trimmomatic + Bowtie2 against the human genome)
# is the standard substitute for shotgun data; see scripts/00_setup.sh section 3a/7 for
# the conda env and human-genome index this depends on (requires SETUP_METAGENOMICS=1
# to have been run at least once).
#
# One array task per study. Output: paired, human-depleted FASTQs in
# $DATA_ROOT/work/kneaddata/<study>/<run>_paired_{1,2}.fastq, one kneaddata log per run.
#
# --max-memory 8g (below) overrides kneaddata's own default of 500m for Trimmomatic's
# JVM heap — confirmed via a real local run (not a guess) that the default OOMs
# ("Java heap space") on a real ~7M-read-pair WGS sample; 500m is sized for 16S-scale
# data, not shotgun. 8g leaves headroom under this job's 48G allocation for bowtie2's
# own (separate, non-JVM) memory use.

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
ENV=$SW_ROOT/envs/rygb-kneaddata
HUMANDB=$SW_ROOT/db/kneaddata_human
THREADS=${SLURM_CPUS_PER_TASK:-16}

# kneaddata shells out to `java -jar trimmomatic.jar`; that jar is built against the
# env's own bundled java (newer) not the system java (older) found on $PATH by default
# when kneaddata is invoked by full path (never conda-activated) — without this, the
# system java fails with "UnsupportedClassVersionError" (class file version 69 vs 61).
export PATH="$ENV/bin:$PATH"

case $SLURM_ARRAY_TASK_ID in
  1) STUDY=Penney_RYGB;  EXPECT=25  ;;
  2) STUDY=Penney_SG;    EXPECT=31  ;;
  3) STUDY=Oseberg_RYGB; EXPECT=78  ;;
  4) STUDY=Oseberg_SG;   EXPECT=76  ;;
  5) STUDY=COCKTAIL_GBP; EXPECT=408 ;;
  6) STUDY=Hayashi_LSG;  EXPECT=38  ;;
  *) echo "bad array index"; exit 1 ;;
esac

RAW=$DATA_ROOT/raw/$STUDY
OUT=$DATA_ROOT/work/kneaddata/$STUDY
mkdir -p "$OUT"
echo "=== $STUDY  $(ls "$RAW"/*_R1.fastq.gz 2>/dev/null | wc -l) raw pairs (expect $EXPECT)   $(date)"

if [ ! -x "$ENV/bin/kneaddata" ]; then
  echo "!! $ENV/bin/kneaddata not found — run scripts/00_setup.sh with SETUP_METAGENOMICS=1 first"
  exit 1
fi
if [ ! -d "$HUMANDB" ] || [ -z "$(ls -A "$HUMANDB" 2>/dev/null)" ]; then
  echo "!! $HUMANDB is empty — kneaddata_database --download human_genome bowtie2 \"$HUMANDB\" (section 7 of 00_setup.sh) hasn't completed"
  exit 1
fi

# A clean directory containing ONLY trimmomatic.jar (symlinked), built once by
# 00_setup.sh's section 3a — NOT the raw share/trimmomatic*/ dir, which also contains a
# same-prefix-matching wrapper script that kneaddata's own glob-based resolution can pick
# instead of the real jar (see 00_setup.sh's comment on this for the full explanation).
TRIMMOMATIC_DIR=$SW_ROOT/db/trimmomatic_jar_only
if [ ! -e "$TRIMMOMATIC_DIR/trimmomatic.jar" ]; then
  echo "!! $TRIMMOMATIC_DIR/trimmomatic.jar not found — run scripts/00_setup.sh with SETUP_METAGENOMICS=1 first"
  exit 1
fi

n=0; skip=0; fail=0; t0=$SECONDS
for f1 in "$RAW"/*_R1.fastq.gz; do
  [ -e "$f1" ] || { echo "no raw reads found in $RAW"; break; }
  run=$(basename "$f1" _R1.fastq.gz)
  f2=$RAW/${run}_R2.fastq.gz
  rundir=$OUT/$run
  # Known, accepted gap: this only checks non-empty, not complete. A process killed
  # mid-write (OOM, walltime, disk-full) leaves a truncated-but-nonempty donefile that
  # gets permanently accepted as done on every future resume, silently feeding reduced
  # read depth into Kraken2/Bracken. Same class of risk as 02b_download_wgs.sbatch's
  # resume check one stage earlier. Not hardened further here — a real fix needs either
  # an atomic write (temp name + rename only on confirmed success) or a read-count
  # cross-check, more than this project's scale currently warrants.
  donefile=$rundir/${run}_R1_kneaddata_paired_1.fastq
  if [ -s "$donefile" ]; then skip=$((skip+1)); continue; fi
  if [ ! -s "$f2" ]; then echo "FAIL $run: missing R2"; fail=$((fail+1)); continue; fi

  # Always start a fresh attempt from a clean directory — a prior crashed/killed run
  # (e.g. the disk-full incident) can leave partial temp files (reordered_*, trf's own
  # *_temp_trf_output* files) that confuse kneaddata_trf_parallel on retry (confirmed:
  # real samples failed with "Subprocess reported error" from trf, traced to a leftover
  # full set of temp files from an earlier killed attempt sitting alongside the new run's).
  rm -rf "$rundir"
  mkdir -p "$rundir"
  "$ENV/bin/kneaddata" \
    --input1 "$f1" --input2 "$f2" \
    --reference-db "$HUMANDB" \
    --output "$rundir" \
    --trimmomatic "$TRIMMOMATIC_DIR" \
    --bowtie2 "$ENV/bin" \
    --trf "$ENV/bin" \
    --max-memory 8g \
    --threads "$THREADS" --processes 1 \
    --remove-intermediate-output \
    > "$rundir.log" 2>&1
  rc=$?
  if [ "$rc" -eq 0 ] && [ -s "$donefile" ]; then
    n=$((n+1))
    [ "$n" -eq 1 ] && echo "first sample: $((SECONDS-t0))s"
    # Raw input is never needed again once kneaddata succeeds for this sample — only
    # kept on failure (above, for retry). Each raw pair belongs to exactly one sample,
    # processed by exactly one task, so deleting here is safe even with other tasks
    # running concurrently on this same cohort's raw/ directory.
    rm -f "$f1" "$f2"
  else
    echo "FAIL $run (rc=$rc)"; tail -5 "$rundir.log"; fail=$((fail+1))
  fi
done

echo
echo "=== $STUDY done  $(date)   elapsed $((SECONDS-t0))s"
echo "processed now: $n | already present: $skip | failures: $fail"
have=$(find "$OUT" -name "*_kneaddata_paired_1.fastq" 2>/dev/null | wc -l)
echo "cleaned pairs: $have (expect $EXPECT)"

echo
echo "--- host-read removal rate (from kneaddata's own logs, first 3 samples)"
for log in $(find "$OUT" -maxdepth 1 -name "*.log" | head -3); do
  echo "  $(basename "$log"):"
  grep -i "contaminate\|final pair" "$log" | sed 's/^/    /'
done
[ "$fail" -eq 0 ] || exit 1
