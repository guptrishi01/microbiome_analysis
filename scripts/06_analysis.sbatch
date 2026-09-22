#!/bin/bash
#SBATCH --job-name=rygb_analysis
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/%u/rygb/logs/analysis_%j.out
#SBATCH --error=/scratch/%u/rygb/logs/analysis_%j.out
#
# Stage 6: the repo's R modules, in the authors' own Singularity images.
#
# The scripts locate their inputs through BioLockJ's directory convention:
#   pipeRoot   = dirname(dirname(getwd()))     <- $PIPE
#   moduleDir  = dirname(getwd())              <- $PIPE/NN_ModuleName
#   and they find upstream modules with dir(pipeRoot, "<ModuleName>")
# So each module gets $PIPE/NN_Name/{script,output,resources} and is run from script/.
#
# Phase 1 = the 16S half: Figures 1, 2, 4 and the diversity/SV supplements.
# Modules needing metagenomic input (OppPathogens, *Metagenomics, *pathway*) are listed
# but skipped until phase 2.

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
REPO=$DATA_ROOT/RYGB_IntegratedAnalysis2020
RSCRIPTS=$REPO/Analysis/RScripts
OUR_RSCRIPTS=$SW_ROOT/analysis/RScripts
PIPE=$DATA_ROOT/pipeline
IMG=$SW_ROOT/images
RUN_PHASE2=${RUN_PHASE2:-0}
# PHASE2_COHORTS=1 adds the 4 phase-2 cohorts (Cholecystectomy, IleocecalResection, Ileostomy,
# SDT) and switches combineCountTables.R/compareStudies_16S.R/compareStudies_SV.R/Heatmap.R to
# this repo's generalized versions (analysis/RScripts/), which cover all 8 studies. Default (0):
# every phase-2 module is skipped and the 4 shared scripts stay the ORIGINAL repo versions,
# copied unmodified from $RSCRIPTS exactly as before — the new code paths below are not run at
# all, not just "verified equivalent," when this flag is off. See CLAUDE.md's "Never assume"
# and "Adding a new cohort" sections: this repo's combineCountTables.R/compareStudies_16S.R/
# Heatmap.R/compareStudies_SV.R were each verified byte-identical to the real phase-1 output
# before being wired in here.
PHASE2_COHORTS=${PHASE2_COHORTS:-0}

BIND="--bind $DATA_ROOT --bind $SW_ROOT"

# name | script path under Analysis/RScripts | image | param | phase
MODULES=(
  "AfsharTaxaClass|16S_TaxaClassification/Afshar_16S_TaxaClassification.R|rbase||1"
  "AssalTaxaClass|16S_TaxaClassification/Assal_16S_TaxaClassification.R|rbase||1"
  "BSTaxaClass|16S_TaxaClassification/BS_16S_TaxaClassification.R|rbase||1"
  "IlhanTaxaClass|16S_TaxaClassification/Ilhan_16S_TaxaClassification.R|rbase||1"
  "CholecystectomyTaxaClass|16S_TaxaClassification/Cholecystectomy_16S_TaxaClassification.R|rbase||2c"
  "IleocecalResectionTaxaClass|16S_TaxaClassification/IleocecalResection_16S_TaxaClassification.R|rbase||2c"
  "IleostomyTaxaClass|16S_TaxaClassification/Ileostomy_16S_TaxaClassification.R|rbase||2c"
  "SDTTaxaClass|16S_TaxaClassification/SDT_16S_TaxaClassification.R|rbase||2c"
  "CombineCountTableslog10|combineCountTables.R|nlme|log10|1"
  "CombineCountTablesrelab|combineCountTables.R|nlme|relab|1"
  "PCO|pco.R|r-vegan||1"
  "AfsharDADA|DADA2_16S/Afshar_DADA2_16S.R|rbase||1"
  "AssalDADA|DADA2_16S/Assal_DADA2_16S.R|rbase||1"
  "BSDADA|DADA2_16S/BS_DADA2_16S.R|rbase||1"
  "IlhanDADA|DADA2_16S/Ilhan_DADA2_16S.R|rbase||1"
  "CholecystectomyDADA|DADA2_16S/Cholecystectomy_DADA2_16S.R|rbase||2c"
  "IleocecalResectionDADA|DADA2_16S/IleocecalResection_DADA2_16S.R|rbase||2c"
  "IleostomyDADA|DADA2_16S/Ileostomy_DADA2_16S.R|rbase||2c"
  "SDTDADA|DADA2_16S/SDT_DADA2_16S.R|rbase||2c"
  "CompareStudies16S|compareStudies/compareStudies_16S.R|comparestudies||1"
  "Heatmap|Heatmap.R|complexheatmap||1"
  "AfsharKraken|kraken2_16S/Afshar_kraken2_16S.R|nlme||1"
  "AssalKraken|kraken2_16S/Assal_kraken2_16S.R|nlme||1"
  "BSKraken|kraken2_16S/BS_kraken2_16S.R|nlme||1"
  "IlhanKraken|kraken2_16S/Ilhan_kraken2_16S.R|nlme||1"
  "CholecystectomyKraken|kraken2_16S/Cholecystectomy_kraken2_16S.R|nlme||2c"
  "IleocecalResectionKraken|kraken2_16S/IleocecalResection_kraken2_16S.R|nlme||2c"
  "IleostomyKraken|kraken2_16S/Ileostomy_kraken2_16S.R|nlme||2c"
  "SDTKraken|kraken2_16S/SDT_kraken2_16S.R|nlme||2c"
  "BSMetagenomics|kraken2_Metagenomics/BS_kraken2_Metagenomics.R|nlme||2"
  "PallejaMetagenomics|kraken2_Metagenomics/Palleja_kraken2_Metagenomics.R|nlme||2"
  "OppPathogens|OpportunisticPathogens.R|pathogens||2"
  "CompareStudiesMetagenomics|compareStudies/compareStudies_Metagenomics.R|comparestudies||2"
  "BSpathway|pathway/BS_pathway.R|nlme||2"
  "Pallejapathway|pathway/Palleja_pathway.R|nlme||2"
  "CompareStudiesPathways|compareStudies/compareStudies_pathways.R|comparestudies||2"
  "Diversity|Diversity.R|r-vegan||1"
  "CompareStudiesSV|compareStudies/compareStudies_SV.R|comparestudies||1"
  "SVAnalysis|SequenceVariantAnalysis.R|comparestudies||1"
  "TrainModellasso|train_model.R|siamcat|lasso|1"
  "FigPerformlasso|Figure_performance.R|fig-perform|lasso|1"
)

echo "=== stage 6  $(date)"
echo "pipeline root: $PIPE"

# ---------------------------------------------------------------- scaffold
mkdir -p "$PIPE"
# the scripts read paste0(pipeRoot,"/input/...") — mirror our regenerated tables there
rm -rf "$PIPE/input"
cp -r "$DATA_ROOT/input" "$PIPE/input"
echo "input tree staged: $(find "$PIPE/input" -type f | wc -l) files"

rm -rf "$PIPE"/[0-9][0-9]_*          # no stale outputs from a previous run
i=0
for entry in "${MODULES[@]}"; do
  IFS='|' read -r name script image param phase <<< "$entry"
  # Skip phase-2-cohort entries BEFORE consuming a number, unlike the pre-existing phase="2"
  # check below (which deliberately keeps its old behavior — see the run loop for why).
  # With PHASE2_COHORTS unset/0 (default), this means every phase-1 module keeps EXACTLY the
  # same numeric prefix already published in results/phase1/pipeline/ — zero directory-naming
  # drift, not just zero content drift.
  [ "$phase" = "2c" ] && [ "$PHASE2_COHORTS" != "1" ] && continue
  num=$(printf "%02d" $i); i=$((i+1))
  [ "$phase" = "2" ] && [ "$RUN_PHASE2" != "1" ] && continue
  d=$PIPE/${num}_${name}
  mkdir -p "$d/script" "$d/output" "$d/resources" "$d/temp"
  if [ "$phase" = "2c" ]; then
    cp "$OUR_RSCRIPTS/$script" "$d/script/"
  else
    cp "$RSCRIPTS/$script" "$d/script/"
  fi
  cp "$RSCRIPTS/functions.R" "$d/resources/"
done
echo "modules scaffolded: $(ls -d "$PIPE"/[0-9][0-9]_* | wc -l)"

# When phase-2 cohorts are included, swap in this repo's generalized versions of the 4 shared
# downstream scripts that hardcode study names in the original (combineCountTables.R,
# compareStudies_16S.R, Heatmap.R, compareStudies_SV.R) — each verified byte-identical to the
# real phase-1 output before being wired in here (see CLAUDE.md's "Never assume" section).
# With PHASE2_COHORTS unset/0 (default), this block does not run — those modules keep the
# original $RSCRIPTS copy from the loop above, untouched.
if [ "$PHASE2_COHORTS" = "1" ]; then
  echo "PHASE2_COHORTS=1: swapping in generalized combineCountTables.R/compareStudies_16S.R/Heatmap.R/compareStudies_SV.R"
  for d in "$PIPE"/[0-9][0-9]_CombineCountTableslog10 "$PIPE"/[0-9][0-9]_CombineCountTablesrelab; do
    [ -d "$d" ] && cp "$OUR_RSCRIPTS/combineCountTables.R" "$d/script/"
  done
  for d in "$PIPE"/[0-9][0-9]_CompareStudies16S; do
    [ -d "$d" ] && cp "$OUR_RSCRIPTS/compareStudies/compareStudies_16S.R" "$d/script/"
  done
  for d in "$PIPE"/[0-9][0-9]_Heatmap; do
    [ -d "$d" ] && cp "$OUR_RSCRIPTS/Heatmap.R" "$d/script/"
  done
  for d in "$PIPE"/[0-9][0-9]_CompareStudiesSV; do
    [ -d "$d" ] && cp "$OUR_RSCRIPTS/compareStudies/compareStudies_SV.R" "$d/script/"
  done
fi

# ---------------------------------------------------------------- convergence guard
# The authors wrapped lme() in tryCatch in BS_kraken2_16S.R ("Streptomycetales and
# Streptomycetaceae produce singular convergence error") but left the other scripts bare —
# with their database, only BS tripped it. The standard Kraken2 index puts taxa in the other
# studies that also fail to converge, so apply their pattern to every copied script.
# Patches the pipeline copies only; the reference clone is never touched.
module load python/3.12.3 2>/dev/null
python3 - "$PIPE" <<'PYEOF'
import re, sys, glob, os
patched = []
for path in glob.glob(os.path.join(sys.argv[1], '[0-9][0-9]_*', 'script', '*.R')):
    src = open(path).read()
    m = re.search(r'^([ \t]*)fit<-anova\((lme\(.*\))\)[ \t]*$', src, re.M)
    if not m:
        continue                      # already guarded (BS) or no model loop
    ind, model = m.group(1), m.group(2)
    new = (f'{ind}fit<-tryCatch({{anova({model})}},'
           f'error=function(e){{cat("LME ERROR:",conditionMessage(e),"\\n");return(NA)}})\n'
           f'{ind}if(is.na(fit)[1]){{ cat("  skipping non-converging column",i,"\\n"); next }}')
    open(path, 'w').write(src[:m.start()] + new + src[m.end():])
    patched.append(os.path.basename(path))
print("convergence guard applied to:", ', '.join(sorted(patched)) if patched else 'nothing')
PYEOF

# ---------------------------------------------------------------- run
# Independent per-cohort modules (name ends in TaxaClass/DADA/Kraken) are batched into groups
# of consecutive array entries sharing the same suffix and run CONCURRENTLY (background jobs,
# capped at SLURM_CPUS_PER_TASK), instead of one at a time. This is safe because each one only
# ever touches its own cohort's input/output paths — confirmed by reading every one of them,
# phase-1 and phase-2 alike — so there is no shared mutable state between concurrent runs.
# Every other module (CombineCountTables, PCO, CompareStudies*, Heatmap, Diversity,
# SVAnalysis, TrainModellasso, FigPerformlasso) still runs alone, exactly as before — an
# untagged module always ends its "batch" at size 1, which is bitwise the same command as the
# old foreground invocation, just plumbed through the same background+exit-code-file mechanism
# for code-path uniformity.
#
# Why this exists: with 8 cpus-per-task requested but each MixedLinearModel R script fitting
# one taxon at a time in a plain sequential for-loop (no parallel/foreach/mclapply), the
# original one-module-at-a-time loop left 7 of 8 requested cores idle for nearly the whole job.
# Batching independent per-cohort modules is the safe fix — it doesn't touch a single line of
# R (the untouchable-clone rule and the "don't restructure someone else's stats code lightly"
# principle both argue against parallelizing *inside* those scripts instead).
#
# Verified, not assumed: see CLAUDE.md's "Never assume" section — this was run against the
# real BS/Assal/Ilhan/Afshar TaxaClass inputs already on disk and diffed byte-for-byte against
# the published Genus_norm_table_*.txt outputs before being trusted here.
MAXPAR=${SLURM_CPUS_PER_TASK:-8}
i=0; ran=0; failed=0
n=${#MODULES[@]}
idx=0
while [ $idx -lt $n ]; do
  batch_nums=(); batch_names=(); batch_dirs=(); wave_tag=""

  while [ $idx -lt $n ] && [ ${#batch_dirs[@]} -lt $MAXPAR ]; do
    entry="${MODULES[$idx]}"
    IFS='|' read -r name script image param phase <<< "$entry"

    # Same skip-before-numbering order as the scaffold loop above — must match exactly, or a
    # module's number here would disagree with the directory the scaffold loop actually created.
    if [ "$phase" = "2c" ] && [ "$PHASE2_COHORTS" != "1" ]; then
      echo "-- skip ${name} (phase 2 cohort: set PHASE2_COHORTS=1 to include)"
      idx=$((idx+1)); continue
    fi

    tag=""
    case "$name" in
      *TaxaClass) tag="TaxaClass" ;;
      *DADA) tag="DADA" ;;
      *Kraken) tag="Kraken" ;;
    esac
    # A later entry with a different tag ends the batch here without consuming it — the next
    # outer-loop iteration picks it up fresh.
    [ ${#batch_dirs[@]} -gt 0 ] && [ "$tag" != "$wave_tag" ] && break
    wave_tag="$tag"

    num=$(printf "%02d" $i); i=$((i+1))
    if [ "$phase" = "2" ] && [ "$RUN_PHASE2" != "1" ]; then
      echo "-- skip ${num}_${name} (phase 2: needs metagenomic input)"
      idx=$((idx+1)); continue
    fi

    d=$PIPE/${num}_${name}
    batch_nums+=("$num"); batch_names+=("$name"); batch_dirs+=("$d")
    base=$(basename "$script")
    echo
    echo "======================================== ${num}_${name}  [$image] $(date +%H:%M:%S)$( [ -n "$tag" ] && echo '  (parallel batch)' )"
    { ( cd "$d/script" && singularity exec $BIND "$IMG/$image.sif" Rscript "./$base" $param ) \
        > "$d/module.log" 2>&1; echo $? > "$d/exit_code"; } &

    idx=$((idx+1))
    [ -z "$tag" ] && break   # untagged modules never batch with anything
  done

  wait
  batch_failed=0
  for k in "${!batch_dirs[@]}"; do
    d="${batch_dirs[$k]}"; nm="${batch_names[$k]}"; nb="${batch_nums[$k]}"
    rc=$(cat "$d/exit_code" 2>/dev/null || echo 1)
    nout=$(find "$d/output" -type f 2>/dev/null | wc -l)
    if [ "$rc" -eq 0 ]; then
      echo "   ok   ${nb}_${nm}  rc=0  outputs=$nout"
      ran=$((ran+1))
    else
      echo "   FAIL ${nb}_${nm}  rc=$rc  outputs=$nout"
      echo "   --- last 25 lines of $d/module.log"
      tail -25 "$d/module.log" | sed 's/^/   /'
      failed=$((failed+1))
      batch_failed=1
    fi
  done
  # dependencies are linear across waves; stop rather than cascade failures into the next one
  [ "$batch_failed" -eq 1 ] && break
done

# ---------------------------------------------------------------- report
echo
echo "=== done  $(date)   modules ok: $ran   failed: $failed"
echo
echo "--- figures produced"
find "$PIPE" -name "*.pdf" -newermt "-1 day" | sort | sed 's/^/   /'
echo
echo "--- key tables vs the repo's published versions"
for f in Genus_PairwiseComparison.txt Genus_countTable_merged_log10.txt metaData_merged.txt; do
  ours=$(find "$PIPE" -name "$f" | head -1)
  theirs=$(ls "$REPO/output/$f" 2>/dev/null || ls "$REPO/output/${f/Comparison/Comaprison}" 2>/dev/null)
  if [ -n "$ours" ] && [ -n "$theirs" ]; then
    printf "   %-38s ours %5s lines | repo %5s lines\n" "$f" "$(wc -l < "$ours")" "$(wc -l < "$theirs")"
  fi
done
echo
echo "--- THE headline number: mean Spearman r across different-study pairs"
echo "    the paper's abstract says 0.41 +/- 0.10; the repo's own table gives 0.410 +/- 0.098"
summarize_r() {
  awk -F'\t' -v tag="$2" 'NR>1{
      gsub(/"/,"",$2); gsub(/"/,"",$5);
      a+=$5; aa+=$5*$5; na++;
      if($2=="Different study"){d+=$5; dd+=$5*$5; nd++}
    }
    END{
      ma=a/na; md=d/nd;
      printf "   %-5s all pairs: r = %.3f +/- %.3f (n=%d) | different-study: r = %.3f +/- %.3f (n=%d)\n",
             tag, ma, sqrt(aa/na-ma*ma), na, md, sqrt(dd/nd-md*md), nd
    }' "$1"
}
ours=$(find "$PIPE" -name "Genus_PairwiseComparison.txt" | head -1)
[ -n "$ours" ] && summarize_r "$ours" "ours"
summarize_r "$REPO/output/Genus_PairwiseComaprison.txt" "repo"
[ "$failed" -eq 0 ] || exit 1
