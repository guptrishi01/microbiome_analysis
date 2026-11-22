#!/bin/bash
#SBATCH --job-name=rygb_k2_wgs
#SBATCH --partition=Orion
#SBATCH --array=1-6
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/%u/rygb/logs/kraken2_wgs_%A_%a.out
#SBATCH --error=/scratch/%u/rygb/logs/kraken2_wgs_%A_%a.out
#
# Stage 4 (WGS track): Kraken2 + Bracken directly on kneaddata-cleaned shotgun reads —
# the SAME k2_standard_20201202 index 00_setup.sh already builds for the 16S track
# (04_kraken2_16s.sbatch), reused as-is. Unlike 16S (where Kraken2 classifies short
# DADA2 ASVs), this runs on full paired-end shotgun reads, which is what a standard
# Kraken2 database is actually built for — no new reference database needed.
#
# Bracken re-estimates abundance at a chosen rank from Kraken2's read-level calls — more
# accurate than reading percentages straight off the Kraken2 report for WGS data. Run at
# every rank so 05b_build_tables_wgs.sbatch can produce the same Phylum..Species table
# set the 16S track and BS-WGS (Fouladi 2021) both already use. Bracken's -r must match a
# k-mer distribution 00_setup.sh's section 6b actually built (150bp for Penney/Oseberg/
# COCKTAIL_GBP, 300bp for Hayashi_LSG's MiSeq 2x300 reads — see the READLEN map below).
#
# Output: one Kraken2 report + one Bracken report per rank per sample, in
# work/kraken_wgs/<study>/<run>.k2report and <run>.<rank>.bracken

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
KDB=$SW_ROOT/db/kraken2
BRACKEN=$SW_ROOT/envs/rygb-kneaddata/bin/bracken
THREADS=${SLURM_CPUS_PER_TASK:-16}
module load kraken2/2.1.3

# bracken's own wrapper script shells out to bare `python $DIR/src/est_abundance.py`
# (not python3, not a full path) — deterministic only if the env's own python (which
# the script's dependencies were resolved against) is first on $PATH, not whatever the
# system/module environment happens to provide.
export PATH="$SW_ROOT/envs/rygb-kneaddata/bin:$PATH"

case $SLURM_ARRAY_TASK_ID in
  1) STUDY=Penney_RYGB;  EXPECT=25  ;;
  2) STUDY=Penney_SG;    EXPECT=31  ;;
  3) STUDY=Oseberg_RYGB; EXPECT=78  ;;
  4) STUDY=Oseberg_SG;   EXPECT=76  ;;
  5) STUDY=COCKTAIL_GBP; EXPECT=408 ;;
  6) STUDY=Hayashi_LSG;  EXPECT=38  ;;
  *) echo "bad array index"; exit 1 ;;
esac

KNEAD=$DATA_ROOT/work/kneaddata/$STUDY
OUT=$DATA_ROOT/work/kraken_wgs/$STUDY
mkdir -p "$OUT"
echo "=== $STUDY  $(find "$KNEAD" -name '*_kneaddata_paired_1.fastq' 2>/dev/null | wc -l) cleaned samples (expect $EXPECT)   $(date)"
echo "db: $KDB  ($(du -sh "$KDB" | cut -f1))"

# Kraken2 reports use single-letter rank codes; map to the same names used downstream
declare -A RANK_NAME=( [P]=Phylum [C]=Class [O]=Order [F]=Family [G]=Genus [S]=Species )

# Per-study read length for Bracken's -r (must match the Bracken k-mer distribution built
# in 00_setup.sh's section 6b — see that file's comment: Bracken's own guidance is ~50bp
# tolerance, and Hayashi (MiSeq 2x300) is 3x over a 150bp distribution, so it needs its own.
declare -A READLEN=( [Penney_RYGB]=150 [Penney_SG]=150 [Oseberg_RYGB]=150 [Oseberg_SG]=150 [COCKTAIL_GBP]=150 [Hayashi_LSG]=300 )
RLEN=${READLEN[$STUDY]}

n=0; skip=0; fail=0; t0=$SECONDS

# Pass 0: backfill any sample that already has a Kraken2 report but is missing one or
# more Bracken rank outputs — a prior Bracken failure for just some ranks. This must run
# as its OWN pass over $OUT, independent of whether that sample's kneaddata directory
# still exists under $KNEAD: the main loop below deletes $d unconditionally once Kraken2
# succeeds (Bracken only ever needs the .k2report it already has, never the original
# reads), so a sample whose Bracken partially failed would otherwise vanish from view on
# every future run (the main loop only iterates over $KNEAD, and that sample's directory
# is already gone) — this was a real gap: previously, only the .k2report gated the
# "already done" skip check, so a Bracken failure would be silently permanent.
for k2out in "$OUT"/*.k2report; do
  [ -s "$k2out" ] || continue
  run=$(basename "$k2out" .k2report)
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
  r1=$d/${run}_R1_kneaddata_paired_1.fastq
  r2=$d/${run}_R1_kneaddata_paired_2.fastq
  [ -s "$r1" ] && [ -s "$r2" ] || continue

  k2out=$OUT/$run.k2report
  if [ -s "$k2out" ]; then skip=$((skip+1)); continue; fi

  if kraken2 --db "$KDB" --threads "$THREADS" --memory-mapping --paired \
             --report "$k2out" --output "$OUT/$run.kraken2.out" \
             "$r1" "$r2" > "$OUT/$run.log" 2>&1; then
    # Bracken per rank, from this one Kraken2 report — no re-classification needed.
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
    # kneaddata's cleaned FASTQs are never needed again once Kraken2 has classified
    # them — the .k2report this just produced is Bracken's own input for every rank
    # (including the Pass-0 backfill above, on a future run), never the reads
    # themselves. $d lives under work/kneaddata/, entirely separate from $OUT
    # (work/kraken_wgs/) where the .k2report/.bracken files this deletes nothing from
    # actually live.
    rm -rf "$d"
  else
    echo "FAIL $run"; tail -3 "$OUT/$run.log"; fail=$((fail+1))
  fi
  rm -f "$OUT/$run.kraken2.out"   # per-read classification, not needed downstream; reports are
done

echo
echo "=== $STUDY done  $(date)   elapsed $((SECONDS-t0))s"
echo "classified now: $n | already present: $skip | failures: $fail"
have=$(ls "$OUT"/*.k2report 2>/dev/null | wc -l)
echo "k2 reports: $have (expect $EXPECT)"

echo
echo "--- classification rates (from the per-sample logs)"
grep -h "sequences classified" "$OUT"/*.log 2>/dev/null \
  | sed 's/.*(\([0-9.]*\)%).*/\1/' | sort -n \
  | awk '{a[NR]=$1} END{if(NR>0) printf "    min %.2f%%  median %.2f%%  max %.2f%%\n", a[1], a[int(NR/2)+1], a[NR]}'
[ "$fail" -eq 0 ] || exit 1
