#!/bin/bash
#SBATCH --job-name=rygb_phase4_real
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/%u/rygb/logs/phase4_real_%j.out
#SBATCH --error=/scratch/%u/rygb/logs/phase4_real_%j.out
#
# Real-data run of the phase-4 WGS R modules (per-cohort MLM + cross-cohort comparison),
# mirroring the earlier fake end-to-end test's structure (scaffold -> convergence-guard
# patch -> run -> report) but against REAL tables built by 05b_build_tables_wgs.sbatch,
# not synthetic data. Deliberately does NOT invoke scripts/06_analysis.sbatch, which
# wipes the shared $DATA_ROOT/pipeline at scaffold time (holding real phase 1-3 results)
# — this scaffolds into an ISOLATED $DATA_ROOT/pipeline_phase4 instead, built from a copy
# of the real $DATA_ROOT/input tree, so the shared pipeline directory is never touched.
#
# Expected outcome, same shape as the fake test (now with real numbers behind it):
#   Penney_RYGBMetagenomics, Penney_SGMetagenomics,
#   Oseberg_RYGBMetagenomics, Oseberg_SGMetagenomics    -> succeed (rc 0)
#   COCKTAIL_GBPMetagenomics, Hayashi_LSGMetagenomics    -> FAIL on purpose (rc != 0),
#                                                            stop() on unresolved Arm
#   CompareStudiesMetagenomicsWGS                        -> succeed (rc 0), comparing
#                                                            whichever of the 4 succeeded,
#                                                            producing the real scatterplot

set -o pipefail
REPO=/users/rgupta25/rygb
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
CLONE_RSCRIPTS=$DATA_ROOT/RYGB_IntegratedAnalysis2020/Analysis/RScripts   # read-only; only read from
OUR_RSCRIPTS=$REPO/analysis/RScripts
PIPE=$DATA_ROOT/pipeline_phase4   # ISOLATED from the real, shared $DATA_ROOT/pipeline
IMG=$SW_ROOT/images
BIND="--bind $DATA_ROOT"

rm -rf "$PIPE"
mkdir -p "$PIPE"
cp -r "$DATA_ROOT/input" "$PIPE/input"   # real tables from 05b; copied, original never touched

MODULES=(
  "Penney_RYGBMetagenomics|kraken2_Metagenomics/Penney_RYGB_kraken2_Metagenomics.R|nlme"
  "Penney_SGMetagenomics|kraken2_Metagenomics/Penney_SG_kraken2_Metagenomics.R|nlme"
  "Oseberg_RYGBMetagenomics|kraken2_Metagenomics/Oseberg_RYGB_kraken2_Metagenomics.R|nlme"
  "Oseberg_SGMetagenomics|kraken2_Metagenomics/Oseberg_SG_kraken2_Metagenomics.R|nlme"
  "COCKTAIL_GBPMetagenomics|kraken2_Metagenomics/COCKTAIL_GBP_kraken2_Metagenomics.R|nlme"
  "Hayashi_LSGMetagenomics|kraken2_Metagenomics/Hayashi_LSG_kraken2_Metagenomics.R|nlme"
  "CompareStudiesMetagenomicsWGS|compareStudies/compareStudies_Metagenomics.R|comparestudies"
)

echo "=== scaffolding real R modules against real input tree ==="
i=0
declare -A DIR_OF
for entry in "${MODULES[@]}"; do
  IFS='|' read -r name script image <<< "$entry"
  num=$(printf "%02d" $i); i=$((i+1))
  d=$PIPE/${num}_${name}
  mkdir -p "$d/script" "$d/output" "$d/resources" "$d/temp"
  cp "$OUR_RSCRIPTS/$script" "$d/script/"
  cp "$CLONE_RSCRIPTS/functions.R" "$d/resources/"
  DIR_OF[$name]=$d
done
echo "modules scaffolded: $(ls -d "$PIPE"/[0-9][0-9]_* | wc -l)"

module load python/3.12.3 2>/dev/null
python3 - "$PIPE" <<'PYEOF'
import re, sys, glob, os
patched = []
for path in glob.glob(os.path.join(sys.argv[1], '[0-9][0-9]_*', 'script', '*.R')):
    src = open(path).read()
    m = re.search(r'^([ \t]*)fit<-anova\((lme\(.*\))\)[ \t]*$', src, re.M)
    if not m:
        continue
    ind, model = m.group(1), m.group(2)
    new = (f'{ind}fit<-tryCatch({{anova({model})}},'
           f'error=function(e){{cat("LME ERROR:",conditionMessage(e),"\\n");return(NA)}})\n'
           f'{ind}if(is.na(fit)[1]){{ cat("  skipping non-converging column",i,"\\n"); next }}')
    open(path, 'w').write(src[:m.start()] + new + src[m.end():])
    patched.append(os.path.basename(os.path.dirname(os.path.dirname(path))))
print("convergence guard applied to:", ', '.join(sorted(patched)) if patched else 'nothing')
PYEOF

echo
echo "=== running real R modules ==="
declare -A RC_OF
for entry in "${MODULES[@]}"; do
  IFS='|' read -r name script image <<< "$entry"
  d=${DIR_OF[$name]}
  base=$(basename "$script")
  echo
  echo "---- $name [$image]"
  ( cd "$d/script" && singularity exec $BIND "$IMG/$image.sif" Rscript "./$base" ) > "$d/module.log" 2>&1
  rc=$?
  RC_OF[$name]=$rc
  nout=$(find "$d/output" -type f 2>/dev/null | wc -l)
  echo "  rc=$rc  outputs=$nout"
  tail -8 "$d/module.log" | sed 's/^/  | /'
done

echo
echo "=== result vs. expectation ==="
declare -A EXPECT=( [Penney_RYGBMetagenomics]=0 [Penney_SGMetagenomics]=0 \
                    [Oseberg_RYGBMetagenomics]=0 [Oseberg_SGMetagenomics]=0 \
                    [COCKTAIL_GBPMetagenomics]=nonzero [Hayashi_LSGMetagenomics]=nonzero \
                    [CompareStudiesMetagenomicsWGS]=0 )
overall_ok=1
for name in "${!EXPECT[@]}"; do
  rc=${RC_OF[$name]}
  exp=${EXPECT[$name]}
  if [ "$exp" = "0" ]; then
    if [ "$rc" -eq 0 ]; then st="PASS"; else st="FAIL"; overall_ok=0; fi
  else
    if [ "$rc" -ne 0 ]; then st="PASS (expected stop())"; else st="FAIL (should have stopped, didn't)"; overall_ok=0; fi
  fi
  printf "  %-32s rc=%-4s expect=%-8s %s\n" "$name" "$rc" "$exp" "$st"
done
echo
if [ "$overall_ok" -eq 1 ]; then
  echo "ALL MODULES MATCHED EXPECTATIONS."
else
  echo "!! MISMATCH — see above. Full logs under $PIPE/*/module.log"
fi

echo
echo "=== CompareStudiesMetagenomicsWGS output (the scatterplot + friends) ==="
find "$PIPE"/*CompareStudiesMetagenomicsWGS/output -type f 2>/dev/null | sed 's/^/  /'
