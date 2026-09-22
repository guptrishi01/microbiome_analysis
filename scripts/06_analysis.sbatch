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
PIPE=$DATA_ROOT/pipeline
IMG=$SW_ROOT/images
RUN_PHASE2=${RUN_PHASE2:-0}

BIND="--bind $DATA_ROOT --bind $SW_ROOT"

# name | script path under Analysis/RScripts | image | param | phase
MODULES=(
  "AfsharTaxaClass|16S_TaxaClassification/Afshar_16S_TaxaClassification.R|rbase||1"
  "AssalTaxaClass|16S_TaxaClassification/Assal_16S_TaxaClassification.R|rbase||1"
  "BSTaxaClass|16S_TaxaClassification/BS_16S_TaxaClassification.R|rbase||1"
  "IlhanTaxaClass|16S_TaxaClassification/Ilhan_16S_TaxaClassification.R|rbase||1"
  "CombineCountTableslog10|combineCountTables.R|nlme|log10|1"
  "CombineCountTablesrelab|combineCountTables.R|nlme|relab|1"
  "PCO|pco.R|r-vegan||1"
  "AfsharDADA|DADA2_16S/Afshar_DADA2_16S.R|rbase||1"
  "AssalDADA|DADA2_16S/Assal_DADA2_16S.R|rbase||1"
  "BSDADA|DADA2_16S/BS_DADA2_16S.R|rbase||1"
  "IlhanDADA|DADA2_16S/Ilhan_DADA2_16S.R|rbase||1"
  "CompareStudies16S|compareStudies/compareStudies_16S.R|comparestudies||1"
  "Heatmap|Heatmap.R|complexheatmap||1"
  "AfsharKraken|kraken2_16S/Afshar_kraken2_16S.R|nlme||1"
  "AssalKraken|kraken2_16S/Assal_kraken2_16S.R|nlme||1"
  "BSKraken|kraken2_16S/BS_kraken2_16S.R|nlme||1"
  "IlhanKraken|kraken2_16S/Ilhan_kraken2_16S.R|nlme||1"
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
  num=$(printf "%02d" $i); i=$((i+1))
  [ "$phase" = "2" ] && [ "$RUN_PHASE2" != "1" ] && continue
  d=$PIPE/${num}_${name}
  mkdir -p "$d/script" "$d/output" "$d/resources" "$d/temp"
  cp "$RSCRIPTS/$script" "$d/script/"
  cp "$RSCRIPTS/functions.R" "$d/resources/"
done
echo "modules scaffolded: $(ls -d "$PIPE"/[0-9][0-9]_* | wc -l)"

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
i=0; ran=0; failed=0
for entry in "${MODULES[@]}"; do
  IFS='|' read -r name script image param phase <<< "$entry"
  num=$(printf "%02d" $i); i=$((i+1))
  if [ "$phase" = "2" ] && [ "$RUN_PHASE2" != "1" ]; then
    echo "-- skip ${num}_${name} (phase 2: needs metagenomic input)"
    continue
  fi
  d=$PIPE/${num}_${name}
  base=$(basename "$script")
  echo
  echo "======================================== ${num}_${name}  [$image] $(date +%H:%M:%S)"
  ( cd "$d/script" && \
    singularity exec $BIND "$IMG/$image.sif" Rscript "./$base" $param ) \
      > "$d/module.log" 2>&1
  rc=$?
  nout=$(find "$d/output" -type f 2>/dev/null | wc -l)
  if [ $rc -eq 0 ]; then
    echo "   ok   rc=0  outputs=$nout"
    ran=$((ran+1))
  else
    echo "   FAIL rc=$rc  outputs=$nout"
    echo "   --- last 25 lines of $d/module.log"
    tail -25 "$d/module.log" | sed 's/^/   /'
    failed=$((failed+1))
    # dependencies are linear here; stop rather than cascade failures
    break
  fi
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
