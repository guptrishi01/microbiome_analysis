#!/bin/bash
#SBATCH --job-name=rygb_kneaddata_oseberg_par
#SBATCH --partition=Orion
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/%u/rygb/logs/kneaddata_oseberg_par_%A_%a.out
#SBATCH --error=/scratch/%u/rygb/logs/kneaddata_oseberg_par_%A_%a.out
#
# One-off parallelized kneaddata for Oseberg_RYGB/Oseberg_SG specifically (76-78 samples
# each). 03b_kneaddata.sbatch runs one cohort per array task, serially looping over every
# sample — fine for Penney_RYGB/SG (25/31 samples, ~7-9h at the measured ~9min/sample),
# but ~21-22h for Oseberg at that scale, too slow for a 12h window. This splits EACH of
# Oseberg_RYGB/SG into 6 parallel chunks (task 1-6 = Oseberg_RYGB chunks 1-6, task 7-12 =
# Oseberg_SG chunks 1-6) by position modulo over the sorted run list — same partitioning
# scheme already used by 02c_resume_cocktail_download.sbatch.
#
# The actual kneaddata invocation (flags, env, fixes) is copied byte-for-byte from
# 03b_kneaddata.sbatch — same --trf/--trimmomatic/--max-memory 8g fixes, verified against
# real data before this was written. If 03b's invocation ever changes, mirror the change
# here too (and vice versa) — there is no single source of truth between them by design,
# this is a one-off, not a permanent pipeline stage; delete once Oseberg_RYGB/SG both
# reach their full expected pair counts.

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
ENV=$SW_ROOT/envs/rygb-kneaddata
HUMANDB=$SW_ROOT/db/kneaddata_human
THREADS=${SLURM_CPUS_PER_TASK:-16}
N_CHUNKS=6

export PATH="$ENV/bin:$PATH"

if [ "$SLURM_ARRAY_TASK_ID" -le "$N_CHUNKS" ]; then
  STUDY=Oseberg_RYGB; CHUNK=$SLURM_ARRAY_TASK_ID; EXPECT=78
else
  STUDY=Oseberg_SG;   CHUNK=$((SLURM_ARRAY_TASK_ID - N_CHUNKS)); EXPECT=76
fi

RAW=$DATA_ROOT/raw/$STUDY
OUT=$DATA_ROOT/work/kneaddata/$STUDY
mkdir -p "$OUT"
echo "=== $STUDY chunk $CHUNK/$N_CHUNKS  $(ls "$RAW"/*_R1.fastq.gz 2>/dev/null | wc -l) raw pairs total (expect $EXPECT)   $(date)"

if [ ! -x "$ENV/bin/kneaddata" ]; then
  echo "!! $ENV/bin/kneaddata not found — run scripts/00_setup.sh with SETUP_METAGENOMICS=1 first"
  exit 1
fi
if [ ! -d "$HUMANDB" ] || [ -z "$(ls -A "$HUMANDB" 2>/dev/null)" ]; then
  echo "!! $HUMANDB is empty — kneaddata_database --download human_genome bowtie2 \"$HUMANDB\" (section 7 of 00_setup.sh) hasn't completed"
  exit 1
fi

TRIMMOMATIC_DIR=$SW_ROOT/db/trimmomatic_jar_only
if [ ! -e "$TRIMMOMATIC_DIR/trimmomatic.jar" ]; then
  echo "!! $TRIMMOMATIC_DIR/trimmomatic.jar not found — run scripts/00_setup.sh with SETUP_METAGENOMICS=1 first"
  exit 1
fi

n=0; skip=0; fail=0; t0=$SECONDS
for f1 in "$RAW"/*_R1.fastq.gz; do
  # A missing $f1 here means either (a) the glob matched nothing at all (nullglob is off,
  # so bash left the literal unexpanded pattern as the one "filename"), or (b) a sibling
  # chunk already deleted its own assigned file between when this task's glob captured
  # the list and when this iteration reached it. Both cases must just skip to the next
  # entry, never break the loop — breaking (the original bug) abandoned every file still
  # further down this task's own captured list, including ones legitimately assigned to
  # THIS chunk (confirmed via a real run, 27899859 task 2: it hit this on a sibling's
  # already-deleted file and silently dropped its own still-unprocessed ERR10905880).
  [ -e "$f1" ] || continue
  run=$(basename "$f1" _R1.fastq.gz)
  # Chunk assignment must be a function of the run's own name, not its position in a
  # directory listing: all N_CHUNKS tasks for this study glob/delete from this same RAW
  # dir concurrently, so a position-based index (idx % N_CHUNKS) drifts between tasks as
  # siblings delete files out from under each other's listing — confirmed via a real run
  # (27885406/27885407) where every one of the 24 array tasks exited 0 but 7 samples were
  # never claimed by any chunk. A hash of the stable run name fixes this.
  h=$(printf '%s' "$run" | cksum | awk '{print $1}')
  [ $(( h % N_CHUNKS + 1 )) -eq "$CHUNK" ] || continue
  f2=$RAW/${run}_R2.fastq.gz
  rundir=$OUT/$run
  donefile=$rundir/${run}_R1_kneaddata_paired_1.fastq
  if [ -s "$donefile" ]; then skip=$((skip+1)); continue; fi
  if [ ! -s "$f2" ]; then echo "FAIL $run: missing R2"; fail=$((fail+1)); continue; fi

  # Always start a fresh attempt from a clean directory — mirrors 03b_kneaddata.sbatch's
  # fix for the same confirmed bug (a prior crashed run's leftover temp files confuse
  # kneaddata_trf_parallel on retry).
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
    # kept on failure (above, for retry). Mirrors 03b_kneaddata.sbatch's cleanup.
    rm -f "$f1" "$f2"
  else
    echo "FAIL $run (rc=$rc)"; tail -5 "$rundir.log"; fail=$((fail+1))
  fi
done

echo
echo "=== $STUDY chunk $CHUNK done  $(date)   elapsed $((SECONDS-t0))s"
echo "processed now: $n | already present: $skip | failures: $fail (this chunk only)"
have=$(find "$OUT" -name "*_kneaddata_paired_1.fastq" 2>/dev/null | wc -l)
echo "total cleaned pairs for $STUDY so far (all chunks): $have (expect $EXPECT)"
[ "$fail" -eq 0 ] || exit 1
