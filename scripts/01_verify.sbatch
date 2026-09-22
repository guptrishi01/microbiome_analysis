#!/bin/bash
#SBATCH --job-name=rygb_verify
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/%u/rygb/logs/verify_%j.out
#SBATCH --error=/scratch/%u/rygb/logs/verify_%j.out

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
ENV16S=$SW_ROOT/envs/rygb-16s
TMP=$DATA_ROOT/work/verify.$SLURM_JOB_ID
mkdir -p "$TMP" "$DATA_ROOT/logs"

sec(){ printf '\n==================== %s\n' "$1"; }
ok(){  printf '   PASS  %s\n' "$1"; }
no(){  printf '   FAIL  %s\n' "$1"; }

sec "node"; hostname; nproc; date

sec "2. DADA2 conda env (binaries called directly, no activation)"
if [ -x "$ENV16S/bin/Rscript" ]; then
  out=$("$ENV16S/bin/Rscript" -e 'suppressMessages(library(dada2)); cat(R.version.string, "| dada2", as.character(packageVersion("dada2")))' 2>&1 | tail -2)
  grep -qi "error" <<<"$out" && { no "dada2"; echo "$out" | sed 's/^/         /'; } || ok "$out"
else
  no "no Rscript at $ENV16S/bin/Rscript"
fi
[ -x "$ENV16S/bin/cutadapt" ] && ok "cutadapt $("$ENV16S/bin/cutadapt" --version 2>&1 | tail -1)" || no "no cutadapt"

sec "3. Kraken2 — load the 47 GB index, classify two real V4 sequences"
module load kraken2/2.1.3
cat > "$TMP/test.fa" <<'FA'
>asv1_expect_Collinsella
TACGTAGGGGGCGAGCGTTATCCGGATTCATTGGGCGTAAAGCGCGCGTAGGCGGCCCGGCAGGCCGGGGGTCGAAGCGGGGGGCTCAACCCCCCGAAGCCCCCGGAACCTCCGCGGCTTGGGTCCGGTAGGGGAGGGTGGAACACCCGGTGTAGCGGTGGAATGCGCAGATATCGGGTGGAACACCGGTGGCGAAGGCG
>asv2_expect_Enterobacteriaceae
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCG
FA
t0=$SECONDS
kraken2 --db "$SW_ROOT/db/kraken2" --threads 8 --use-names \
        --report "$TMP/test.report" --output "$TMP/test.out" "$TMP/test.fa" 2>&1 | sed 's/^/   /'
echo "   index load + classify: $((SECONDS-t0))s"
[ -s "$TMP/test.out" ] && { ok "kraken2 classified"; cut -f1-3 "$TMP/test.out" | sed 's/^/         /'; } || no "no kraken2 output"
echo "   --- mpa-style report (format the table builder must parse)"
kraken2 --db "$SW_ROOT/db/kraken2" --threads 8 --use-mpa-style \
        --report "$TMP/test.mpa" --output /dev/null "$TMP/test.fa" >/dev/null 2>&1
head -10 "$TMP/test.mpa" | sed 's/^/         /'
echo "   --- standard report"
head -8 "$TMP/test.report" | sed 's/^/         /'

sec "4. Network from this compute node"
for u in https://www.ebi.ac.uk https://sra-download.ncbi.nlm.nih.gov https://sra-pub-run-odp.s3.amazonaws.com; do
  printf '   %s %s\n' "$(curl -s -o /dev/null -m 20 -w '%{http_code}' "$u")" "$u"
done
echo "   --- which BS BioProject is the 16S set?"
for p in PRJNA668357 PRJNA668472; do
  tsv=$(curl -s -m 60 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$p&result=read_run&fields=run_accession,sample_alias,sample_title,library_strategy,library_layout&format=tsv")
  n=$(printf '%s\n' "$tsv" | tail -n +2 | grep -c .)
  strat=$(printf '%s\n' "$tsv" | tail -n +2 | cut -f4 | sort -u | tr '\n' ',')
  echo "   $p: $n runs | strategy: $strat"
  printf '%s\n' "$tsv" | sed -n '2,4p' | cut -f1-3 | sed 's/^/        /'
done

sec "5. SRA download path — one small Ilhan run"
module load sra-tools/3.1.0
cd "$TMP" || exit 1
if prefetch --max-size 5G -O "$TMP" SRR11292124 >"$TMP/prefetch.log" 2>&1 && \
   fasterq-dump --split-files -e 8 -O "$TMP" "$TMP/SRR11292124" >"$TMP/fqdump.log" 2>&1; then
  ok "prefetch + fasterq-dump"
  ls -lh "$TMP"/SRR11292124*.fastq | sed 's/^/         /'
  head -2 "$TMP"/SRR11292124_1.fastq | cut -c1-80 | sed 's/^/         /'
else
  no "SRA download failed — downloads move to the DTN partition"
  tail -5 "$TMP/prefetch.log" "$TMP/fqdump.log" 2>/dev/null | sed 's/^/         /'
fi

sec "summary"
df -hP "$DATA_ROOT" | tail -1 | sed 's/^/   /'
echo "   sacct -j $SLURM_JOB_ID --format=JobID,State,Elapsed,MaxRSS,ReqMem"
