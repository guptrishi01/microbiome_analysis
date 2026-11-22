#!/bin/bash
#SBATCH --job-name=rygb_k2_oseberg_par
#SBATCH --partition=Orion
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/%u/rygb/logs/kraken2_oseberg_par_%A_%a.out
#SBATCH --error=/scratch/%u/rygb/logs/kraken2_oseberg_par_%A_%a.out
#
# One-off parallelized Kraken2+Bracken for Oseberg_RYGB/Oseberg_SG — same reasoning as
# 03c_kneaddata_oseberg_parallel.sbatch (04b_kraken2_wgs.sbatch's one-task-per-cohort,
# serial-per-sample loop would take ~10-11h for 76-78 samples; too slow for a 12h window
# alongside the kneaddata stage ahead of it). Splits each study into 6 parallel chunks
# (task 1-6 = Oseberg_RYGB, task 7-12 = Oseberg_SG) by the same position-modulo scheme as
# 03c/02c. Classification logic (flags, env, fixes) copied byte-for-byte from
# 04b_kraken2_wgs.sbatch — no single source of truth between them by design, this is a
# one-off; delete once Oseberg_RYGB/SG both reach their full expected .k2report counts.

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
KDB=$SW_ROOT/db/kraken2
BRACKEN=$SW_ROOT/envs/rygb-kneaddata/bin/bracken
THREADS=${SLURM_CPUS_PER_TASK:-16}
module load kraken2/2.1.3
export PATH="$SW_ROOT/envs/rygb-kneaddata/bin:$PATH"

N_CHUNKS=6
if [ "$SLURM_ARRAY_TASK_ID" -le "$N_CHUNKS" ]; then
  STUDY=Oseberg_RYGB; CHUNK=$SLURM_ARRAY_TASK_ID; EXPECT=78; RLEN=150
else
  STUDY=Oseberg_SG;   CHUNK=$((SLURM_ARRAY_TASK_ID - N_CHUNKS)); EXPECT=76; RLEN=150
fi

KNEAD=$DATA_ROOT/work/kneaddata/$STUDY
OUT=$DATA_ROOT/work/kraken_wgs/$STUDY
mkdir -p "$OUT"
echo "=== $STUDY chunk $CHUNK/$N_CHUNKS  $(find "$KNEAD" -name '*_kneaddata_paired_1.fastq' 2>/dev/null | wc -l) cleaned samples total (expect $EXPECT)   $(date)"
echo "db: $KDB  ($(du -sh "$KDB" | cut -f1))"

declare -A RANK_NAME=( [P]=Phylum [C]=Class [O]=Order [F]=Family [G]=Genus [S]=Species )

n=0; skip=0; fail=0; t0=$SECONDS

# Pass 0: backfill any sample THIS CHUNK owns that already has a Kraken2 report but is
# missing one or more Bracken rank outputs — mirrors 04b_kraken2_wgs.sbatch's same fix
# (a prior Bracken failure for just some ranks used to be permanently invisible, since
# only .k2report gated the "already done" skip and $d gets deleted regardless of
# Bracken's outcome). Hashed the same way as the main loop below so concurrent chunks
# never race on the same k2out file's missing ranks.
for k2out in "$OUT"/*.k2report; do
  [ -s "$k2out" ] || continue
  run=$(basename "$k2out" .k2report)
  h=$(printf '%s' "$run" | cksum | awk '{print $1}')
  [ $(( h % N_CHUNKS + 1 )) -eq "$CHUNK" ] || continue
  bracken_ok=1
  for code in P C O F G S; do
    rank=${RANK_NAME[$code]}
    [ -s "$OUT/$run.$rank.bracken" ] && continue
    "$BRACKEN" -d "$KDB" -i "$k2out" \
      -o "$OUT/$run.$rank.bracken" -w "$OUT/$run.$rank.breport" \
      -r "$RLEN" -l "$code" > "$OUT/$run.$rank.bracken.log" 2>&1 \
      || { echo "  !! bracken backfill failed for $run at rank $rank"; bracken_ok=0; }
  done
  [ "$bracken_ok" -eq 1 ] || fail=$((fail+1))
done

for d in "$KNEAD"/*/; do
  run=$(basename "$d")
  # Hash the stable run name rather than its position in a live directory listing — see
  # 03c_kneaddata_oseberg_parallel.sbatch's comment on this same fix for why: concurrent
  # tasks deleting from this same KNEAD dir make a position-based index drift between
  # tasks, letting samples fall through unclaimed by any chunk (confirmed via a real run).
  h=$(printf '%s' "$run" | cksum | awk '{print $1}')
  [ $(( h % N_CHUNKS + 1 )) -eq "$CHUNK" ] || continue

  r1=$d/${run}_R1_kneaddata_paired_1.fastq
  r2=$d/${run}_R1_kneaddata_paired_2.fastq
  [ -s "$r1" ] && [ -s "$r2" ] || continue

  k2out=$OUT/$run.k2report
  if [ -s "$k2out" ]; then skip=$((skip+1)); continue; fi

  if kraken2 --db "$KDB" --threads "$THREADS" --memory-mapping --paired \
             --report "$k2out" --output "$OUT/$run.kraken2.out" \
             "$r1" "$r2" > "$OUT/$run.log" 2>&1; then
    bracken_ok=1
    for code in P C O F G S; do
      rank=${RANK_NAME[$code]}
      "$BRACKEN" -d "$KDB" -i "$k2out" \
        -o "$OUT/$run.$rank.bracken" -w "$OUT/$run.$rank.breport" \
        -r "$RLEN" -l "$code" > "$OUT/$run.$rank.bracken.log" 2>&1 \
        || { echo "  !! bracken failed for $run at rank $rank"; bracken_ok=0; }
    done
    if [ "$bracken_ok" -eq 1 ]; then
      n=$((n+1))
      [ "$n" -eq 1 ] && echo "first sample (includes index load): $((SECONDS-t0))s"
    else
      fail=$((fail+1))
    fi
    # kneaddata's cleaned FASTQs are never needed again once Kraken2 has classified them
    # — mirrors 04b_kraken2_wgs.sbatch's cleanup. $d (work/kneaddata/) is entirely
    # separate from $OUT (work/kraken_wgs/), so this never touches .k2report/.bracken.
    rm -rf "$d"
  else
    echo "FAIL $run"; tail -3 "$OUT/$run.log"; fail=$((fail+1))
  fi
  rm -f "$OUT/$run.kraken2.out"
done

echo
echo "=== $STUDY chunk $CHUNK done  $(date)   elapsed $((SECONDS-t0))s"
echo "classified now: $n | already present: $skip | failures: $fail (this chunk only)"
have=$(ls "$OUT"/*.k2report 2>/dev/null | wc -l)
echo "total k2 reports for $STUDY so far (all chunks): $have (expect $EXPECT)"
[ "$fail" -eq 0 ] || exit 1
