#!/bin/bash
#SBATCH --job-name=rygb_auroc
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/%u/rygb/logs/auroc_%j.out
#SBATCH --error=/scratch/%u/rygb/logs/auroc_%j.out
#
# Figure 4 as numbers: rebuild the AUROC matrix from our lasso predictions and from the
# repo's published predictions, and print them side by side.

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
REPO=$DATA_ROOT/RYGB_IntegratedAnalysis2020
PIPE=$DATA_ROOT/pipeline

OURS_PRED=$(find "$PIPE" -name "predictions_lasso.tsv" | head -1)
OURS_META=$(find "$PIPE" -path "*CombineCountTablesrelab*" -name "metaData_merged.txt" | head -1)
REPO_PRED=$REPO/output/lasso/predictions_lasso.tsv
REPO_META=$REPO/output/metaData_merged.txt

echo "ours: $OURS_PRED"
echo "      $OURS_META"
echo "repo: $REPO_PRED"

cat > "$DATA_ROOT/work/auroc.R" <<'RSCRIPT'
suppressMessages({library(pROC); library(tidyverse)})
a <- commandArgs(trailingOnly=TRUE)

auroc_matrix <- function(predfile, metafile, label) {
  meta <- read.table(metafile, sep="\t", header=TRUE)
  pred <- read.table(predfile, sep='\t', check.names=FALSE)
  pred$Sample_ID <- rownames(pred)
  df <- inner_join(meta, as_tibble(pred), by='Sample_ID')
  studies <- sort(unique(as.character(meta$Study)))
  m <- matrix(NA, nrow=length(studies), ncol=length(studies)+1,
              dimnames=list(studies, c(studies,'LOSO')))
  for (tr in studies) for (te in studies) {
    d <- df %>% filter(Study == te)
    m[tr, te] <- as.numeric(suppressMessages(
      roc(predictor=d[[tr]], response=d$time, quiet=TRUE)$auc))
  }
  for (te in studies) {
    d <- df %>% filter(Study == te)
    m[te, 'LOSO'] <- as.numeric(suppressMessages(
      roc(predictor=d[['LOSO']], response=d$time, quiet=TRUE)$auc))
  }
  cat("\n===", label, "\n")
  print(round(m, 3))
  diag_auc <- sapply(studies, function(s) m[s, s])
  off <- m[studies, studies]; diag(off) <- NA
  cat(sprintf("  within-study (CV) mean AUROC : %.3f\n", mean(diag_auc)))
  cat(sprintf("  study-to-study transfer mean : %.3f\n", mean(off, na.rm=TRUE)))
  cat(sprintf("  leave-one-study-out mean     : %.3f\n", mean(m[, 'LOSO'], na.rm=TRUE)))
  invisible(m)
}

ours <- auroc_matrix(a[1], a[2], "OURS (regenerated from raw reads)")
repo <- auroc_matrix(a[3], a[4], "REPO (published)")

common_r <- intersect(rownames(ours), rownames(repo))
common_c <- intersect(colnames(ours), colnames(repo))
d <- ours[common_r, common_c] - repo[common_r, common_c]
cat("\n=== difference (ours - repo)\n")
print(round(d, 3))
cat(sprintf("\nmean |difference| = %.3f   max |difference| = %.3f\n",
            mean(abs(d), na.rm=TRUE), max(abs(d), na.rm=TRUE)))
RSCRIPT

singularity exec --bind "$DATA_ROOT" --bind "$SW_ROOT" "$SW_ROOT/images/fig-perform.sif" \
  Rscript "$DATA_ROOT/work/auroc.R" "$OURS_PRED" "$OURS_META" "$REPO_PRED" "$REPO_META"
