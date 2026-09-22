#!/bin/bash
#SBATCH --job-name=rygb_dl
#SBATCH --partition=Orion
#SBATCH --array=1-4
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/%u/rygb/logs/download_%A_%a.out
#SBATCH --error=/scratch/%u/rygb/logs/download_%A_%a.out
#
# Phase 1 downloads: the four 16S studies, forward reads only.
#   task 1 BS16S   PRJNA668472  61 runs   -> renamed to SAMPLE-N (repo's sample IDs)
#   task 2 Assal   PRJNA395660  71 runs   -> kept as SRR accessions (what the repo uses)
#   task 3 Afshar  PRJEB14082   49 runs   -> kept as ERR accessions
#   task 4 Ilhan   PRJNA612042  82 runs   -> kept as SRR accessions
# Resumable: a run whose .fastq.gz already exists is skipped.
# Submit:  sbatch 02_download.sbatch

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
REPO=$DATA_ROOT/RYGB_IntegratedAnalysis2020
module load sra-tools/3.1.0

case $SLURM_ARRAY_TASK_ID in
  1) STUDY=BS16S;  EXPECT=61 ;;
  2) STUDY=Assal;  EXPECT=71 ;;
  3) STUDY=Afshar; EXPECT=49 ;;
  4) STUDY=Ilhan;  EXPECT=82 ;;
  *) echo "bad array index"; exit 1 ;;
esac

OUT=$DATA_ROOT/raw/$STUDY
TMP=$DATA_ROOT/work/dl_$STUDY
mkdir -p "$OUT" "$TMP"
cd "$TMP" || exit 1
echo "=== $STUDY  expecting $EXPECT runs   $(date)"

# ---------------------------------------------------------------- run list
# Two columns: run accession <TAB> output sample name.
LIST=$OUT/runs.tsv
if [ -s "$LIST" ]; then
  echo "run list exists: $LIST"
else
  case $STUDY in
    BS16S)
      # Not in the repo — the repo's metadata carries SAMPLE-N with no accessions.
      # ENA gives alias Sample<N>; normalize to the repo's SAMPLE-<N>.
      curl -s -m 120 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA668472&result=read_run&fields=run_accession,sample_alias&format=tsv" \
        | tail -n +2 \
        | awk -F'\t' '{n=$2; gsub(/[^0-9]/,"",n); if(n!="") print $1"\tSAMPLE-"n}' \
        | sort -t- -k2,2n > "$LIST"
      ;;
    Assal)   # CSV, column 1 = Run
      tail -n +2 "$REPO/input/RYGB_Assal2020/Metadata/SraRunTable_Assal.txt" \
        | cut -d, -f1 | grep . | awk '{print $1"\t"$1}' > "$LIST" ;;
    Afshar)  # TSV, column 1 = Run
      tail -n +2 "$REPO/input/RYGB_Afshar2018/Metadata/metaData.txt" \
        | cut -f1 | grep . | awk '{print $1"\t"$1}' > "$LIST" ;;
    Ilhan)   # TSV whose header is shifted; data column 1 = Run
      tail -n +2 "$REPO/input/RYGB_Ilhan2020/Metadata/metaData_Ilhan.txt" \
        | cut -f1 | grep . | awk '{print $1"\t"$1}' > "$LIST" ;;
  esac
fi

N=$(wc -l < "$LIST")
echo "run list: $N entries (expected $EXPECT)"
[ "$N" -eq "$EXPECT" ] || echo "!! WARNING: count mismatch — check $LIST before trusting downstream output"
head -3 "$LIST"

# ---------------------------------------------------------------- download
fail=0; done_n=0; skip=0
while IFS=$'\t' read -r run sample; do
  target=$OUT/${sample}_R1.fastq.gz
  if [ -s "$target" ]; then skip=$((skip+1)); continue; fi

  if ! prefetch --max-size 100G -O "$TMP" "$run" >"$TMP/$run.prefetch.log" 2>&1; then
    echo "FAIL prefetch $run"; tail -3 "$TMP/$run.prefetch.log"; fail=$((fail+1)); continue
  fi
  if ! fasterq-dump --split-files --threads 8 -O "$TMP" "$TMP/$run" >"$TMP/$run.fqdump.log" 2>&1; then
    echo "FAIL fasterq-dump $run"; tail -3 "$TMP/$run.fqdump.log"; fail=$((fail+1)); rm -rf "$TMP/$run"; continue
  fi

  # forward reads only; single-end submissions land without the _1 suffix
  fwd=$TMP/${run}_1.fastq
  [ -s "$fwd" ] || fwd=$TMP/${run}.fastq
  if [ -s "$fwd" ]; then
    "$SW_ROOT/envs/rygb-16s/bin/pigz" -p 8 -c "$fwd" > "$target"
    done_n=$((done_n+1))
  else
    echo "FAIL no forward reads for $run"; fail=$((fail+1))
  fi
  rm -rf "$TMP/$run" "$TMP/${run}"_*.fastq "$TMP/${run}.fastq"
done < "$LIST"

# ---------------------------------------------------------------- report
echo
echo "=== $STUDY done  $(date)"
echo "downloaded now: $done_n | already present: $skip | failures: $fail"
have=$(ls "$OUT"/*_R1.fastq.gz 2>/dev/null | wc -l)
echo "files in $OUT: $have (expected $EXPECT)"
du -sh "$OUT"
echo "read length check:"
zcat "$(ls "$OUT"/*_R1.fastq.gz | head -1)" | head -2 | tail -1 | awk '{print "   first read: "length($0)" bp: "substr($0,1,50)"..."}'
echo "reads per sample (first 5):"
for f in $(ls "$OUT"/*_R1.fastq.gz | head -5); do
  echo "   $(basename "$f") $(( $(zcat "$f" | wc -l) / 4 )) reads"
done
rm -rf "$TMP"
[ "$fail" -eq 0 ] || exit 1
