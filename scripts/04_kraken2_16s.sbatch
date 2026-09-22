#!/bin/bash
#SBATCH --job-name=rygb_k2
#SBATCH --partition=Orion
#SBATCH --array=1-8
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/%u/rygb/logs/kraken2_%A_%a.out
#SBATCH --error=/scratch/%u/rygb/logs/kraken2_%A_%a.out
#
# Kraken2 on the DADA2-FILTERED forward reads — the same input the authors used
# (their configs point at .../dada2Filtered/bljInput).
# One array task per study; the index is memory-mapped so it loads once and then
# stays in page cache for the rest of the samples.
# Output: one mpa-style report per sample in work/kraken/<study>/<sample>.mpa

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
KDB=$SW_ROOT/db/kraken2
THREADS=${SLURM_CPUS_PER_TASK:-16}
module load kraken2/2.1.3

case $SLURM_ARRAY_TASK_ID in
  1) STUDY=BS16S;              EXPECT=61  ;;
  2) STUDY=Assal;               EXPECT=71  ;;
  3) STUDY=Afshar;              EXPECT=50  ;;
  4) STUDY=Ilhan;               EXPECT=82  ;;
  5) STUDY=Cholecystectomy;     EXPECT=30  ;;
  6) STUDY=IleocecalResection;  EXPECT=189 ;;
  7) STUDY=Ileostomy;           EXPECT=65  ;;
  8) STUDY=SDT;                 EXPECT=45  ;;
  *) echo "bad array index"; exit 1 ;;
esac

FILT=$DATA_ROOT/work/filtered/$STUDY
OUT=$DATA_ROOT/work/kraken/$STUDY
mkdir -p "$OUT"
echo "=== $STUDY  $(ls "$FILT"/*_R1.fastq.gz | wc -l) filtered samples (expect $EXPECT)   $(date)"
echo "db: $KDB  ($(du -sh "$KDB" | cut -f1))"

n=0; skip=0; fail=0; t0=$SECONDS
for f in "$FILT"/*_R1.fastq.gz; do
  s=$(basename "$f" _R1.fastq.gz)
  if [ -s "$OUT/$s.mpa" ]; then skip=$((skip+1)); continue; fi
  if kraken2 --db "$KDB" --threads "$THREADS" --memory-mapping --gzip-compressed \
             --use-mpa-style --report "$OUT/$s.mpa" --output /dev/null "$f" \
             > "$OUT/$s.log" 2>&1; then
    n=$((n+1))
    [ "$n" -eq 1 ] && echo "first sample (includes index load): $((SECONDS-t0))s"
  else
    echo "FAIL $s"; tail -3 "$OUT/$s.log"; fail=$((fail+1))
  fi
done

echo
echo "=== $STUDY done  $(date)   elapsed $((SECONDS-t0))s"
echo "classified now: $n | already present: $skip | failures: $fail"
have=$(ls "$OUT"/*.mpa 2>/dev/null | wc -l)
echo "mpa reports: $have (expect $EXPECT)"

echo
echo "--- classification rates (from the per-sample logs)"
grep -h "sequences classified" "$OUT"/*.log \
  | sed 's/.*(\([0-9.]*\)%).*/\1/' | sort -n \
  | awk '{a[NR]=$1} END{print "    min "a[1]"%  median "a[int(NR/2)+1]"%  max "a[NR]"%"}'

echo "--- a sample report, deepest lines"
f=$(ls "$OUT"/*.mpa | head -1)
echo "    $(basename "$f")"
grep "s__" "$f" | sort -k2 -nr | head -5 | sed 's/^/    /'
echo "--- rank coverage in that sample (how many taxa are named at each rank)"
for r in p c o f g s; do
  printf "    %s__ %s\n" "$r" "$(grep -c "|${r}__[^|]*$" "$f")"
done
[ "$fail" -eq 0 ] || exit 1
