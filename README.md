# Gut Microbiome Signature Correlation Analysis

Fodor Lab project reproducing and extending:

> Fouladi F, et al. **"A microbial signature following bariatric surgery is
> robustly consistent across multiple cohorts."** *Gut Microbes* / *Genome
> Medicine* (2021). Original analysis code:
> [FarnazFouladi/RYGB_IntegratedAnalysis2020](https://github.com/FarnazFouladi/RYGB_IntegratedAnalysis2020)

This README describes what the pipeline is and how to run it. **For actual
findings, numbers, and interpretation, see [results/README.md](results/README.md)
and the per-phase output under `results/`.**

## Objective

1. **Reproduce** the paper's finding — that bariatric surgery (RYGB and
   related procedures) leaves a gut-microbial signature that is consistent
   across independently collected cohorts — by re-running its analysis
   pipeline from raw SRA reads through the same DADA2 / Kraken2 / cross-study
   LASSO steps the authors used, rather than from their pre-computed tables.
2. **Extend** the question beyond bariatric surgery: identify additional
   public 16S datasets from *other* GI-surgery contexts (colectomy,
   gastrectomy, ileostomy, etc.) and test whether a comparable
   surgery-associated microbial signature appears in cohorts the original
   paper never looked at.
3. **Positive-control** the result from (2): identify additional bariatric
   cohorts, split into RYGB-specific and non-RYGB-bariatric sub-phases, to
   test whether the signature is bariatric-specific, RYGB-specific, or
   neither.

Phase 1 (reproduction), Phase 2 (extension to non-bariatric GI surgery), and
Phase 3 (bariatric positive control, sub-phases 3a/3b) are all complete. See
[results/README.md](results/README.md) for what each phase found.

## Repo structure

```
microbiome_analysis/
├── scripts/                  # Pipeline stages, run in order via sbatch (SLURM)
│   ├── 00_setup.sh              # One-time env/software/reference-data setup (idempotent)
│   ├── 01_verify.sbatch         # Verify accessions/metadata for the 4 source cohorts
│   ├── 02_download.sbatch       # Download raw 16S reads from SRA (forward reads only)
│   ├── 03_dada2.sbatch          # DADA2 ASV inference per cohort
│   ├── 04_kraken2_16s.sbatch    # Kraken2 taxonomic classification of DADA2-filtered reads
│   ├── 05_build_tables.sbatch   # Build taxon count tables in the original repo's layout
│   ├── 06_analysis.sbatch       # Run the paper's own R modules (Singularity images)
│   ├── 07_compare_auroc.sbatch  # Compare our LASSO AUROC matrix to the published one
│   ├── 09_find_datasets.sbatch  # Scan ENA for candidate non-bariatric GI-surgery 16S studies
│   └── 10_inspect_candidates.sh # Deep-inspect scan candidates (subject/timepoint/read length)
├── dataset_scan/              # Cached ENA metadata (one TSV per BioProject) + shortlist.tsv
│                               # ranking candidate cohorts for phase 2
├── cohorts/                   # Phase-2/3 per-cohort metadata, built by hand from ENA run info
│   ├── Cholecystectomy/metaData.txt    # phase 2
│   ├── IleocecalResection/metaData.txt # phase 2
│   ├── Ileostomy/metaData.txt          # phase 2
│   ├── SDT/metaData.txt                # phase 2
│   ├── RYGB_SaoPaulo/metaData.txt      # phase 3b (RYGB-specific)
│   ├── Lazaro_GB/metaData.txt          # phase 3b (RYGB-specific)
│   ├── Lazaro_SG/metaData.txt          # phase 3a (non-RYGB bariatric)
│   └── Chaudhari_SG/metaData.txt       # phase 3a (non-RYGB bariatric)
├── analysis/RScripts/         # This repo's own R code, run inside the authors' Singularity
│   │                           # images (never inside their cloned repo)
│   ├── 16S_TaxaClassification/   # One script per phase-2/3 cohort, modeled on the closest
│   ├── DADA2_16S/                 # matching original template (Ilhan/Afshar/Assal/BS —
│   ├── kraken2_16S/                see CLAUDE.md's "Adding a new cohort" section)
│   ├── combineCountTables.R      # Generalized rewrites of the 4 scripts that hardcode all
│   ├── Heatmap.R                  # 4 phase-1 study names in the original repo — verified
│   ├── compareStudies/             # byte-identical to the real phase-1 output before trusting
│   ├── pairwiseHeatmap.R         # Study-timepoint x study-timepoint heatmap (r + BH-adjusted
│   │                               # p per cell) from any run's Genus_PairwiseComparison.txt
│   ├── crossGroupScatter.R       # Re-plots only cross-group pairwise scatterplots (e.g. new
│   │                               # cohorts vs. OG), reusing already-gathered MLM results
│   └── titlePage.R               # One-page section-divider generator for combined reports
├── results/                   # All findings — see results/README.md
│   ├── README.md                 # Findings and interpretation for every phase
│   ├── phase1/                   # Reproduction outputs (figures/, pipeline/, predictions, logs)
│   ├── phase2/                   # Extension outputs, all 8 cohorts
│   └── phase3/                   # Positive-control outputs
│       ├── OG_plus_3a/              # OG + Chaudhari_SG + Lazaro_SG (non-RYGB bariatric) run
│       ├── OG_plus_3b/              # OG + RYGB_SaoPaulo + Lazaro_GB (RYGB-specific) run
│       ├── all12/                   # OG + all phase-2 + all phase-3 cohorts, one combined run
│       ├── report/                  # Combined pairwise-comparison PDF report + r/p CSVs
│       └── *.out                    # Shared download/DADA2/Kraken2/table-building logs
├── setup.log                  # Output of scripts/00_setup.sh (environment provenance)
├── versions.txt                # Tool/container/environment versions used
└── .gitignore                  # Excludes conda_pkgs/, envs/, images/, db/, singularity_cache/
                                 #   (~60GB of downloaded software/reference data — rebuild with
                                 #   scripts/00_setup.sh rather than storing in git)
```

Not tracked in git (built locally by `scripts/00_setup.sh`, `$SW_ROOT`):
`conda_pkgs/`, `envs/` (conda environments), `images/` (the authors'
Singularity/R containers), `db/` (SILVA 132 + Kraken2 index), `singularity_cache/`.

## Workflow (phase 1 — reproduction)

The four cohorts from the original paper — **BS** (the authors' own bariatric
surgery cohort), **Assal2020**, **Afshar2018**, and **Ilhan2020** — are
processed end to end from raw SRA reads rather than starting from the
original repo's precomputed count tables:

1. **`00_setup.sh`** — one-time setup: conda env for DADA2/cutadapt, the
   authors' published Singularity containers (`rbase`, `nlme`, `r-vegan`,
   `comparestudies`, `complexheatmap`, `pathogens`, `siamcat`, `fig-perform`),
   SILVA 132 training sets, a Kraken2 index (the authors' original custom DB
   no longer exists, so a public stand-in of comparable vintage is used),
   and a clone of the original repo for its R modules and sample metadata.
2. **`01_verify.sbatch`** — confirm SRA accessions and metadata for all 4
   cohorts before committing compute to downloads.
3. **`02_download.sbatch`** — download raw forward reads per cohort (SLURM
   array, one task per cohort).
4. **`03_dada2.sbatch`** — DADA2 ASV inference per cohort (same `truncLen`
   convention as the original: 200bp for BS/Assal/Afshar, 150bp for Ilhan
   to match its shorter reads).
5. **`04_kraken2_16s.sbatch`** — Kraken2 classification of the DADA2-filtered
   reads (the same input the authors classified, not raw reads).
6. **`05_build_tables.sbatch`** — convert Kraken2/mpa reports into taxon count
   tables in the original repo's exact file layout and naming convention, so
   the unmodified downstream R modules can consume them.
7. **`06_analysis.sbatch`** — run the original repo's own R analysis modules,
   unmodified, inside the authors' published Singularity images: per-cohort
   taxa classification, combined log10/relative-abundance count tables,
   PCO, cross-study correlation and heatmap figures, diversity, and
   cross-study LASSO training with leave-one-study-out (LOSO) validation.
8. **`07_compare_auroc.sbatch`** — rebuild the LASSO AUROC matrix from our
   regenerated predictions and from the original repo's published
   predictions, and diff them directly.

## Workflow (phase 2 — extending beyond bariatric surgery)

`scripts/09_find_datasets.sbatch` scans ENA study metadata for candidate
non-bariatric GI-surgery 16S cohorts against explicit inclusion criteria
(Illumina amplicon sequencing, ≥150bp reads, human gut samples, repeat
sampling per subject, a recoverable pre-surgical baseline plus ≥1
post-surgical timepoint, GI-surgery context), scoring and ranking results
into `dataset_scan/shortlist.tsv`. `scripts/10_inspect_candidates.sh` then
deep-inspects specific shortlisted accessions — sample label parsing,
subject/timepoint layout, read length — to decide which are actually usable
and with what DADA2 `truncLen`.

### Cohorts used

Four cohorts (three source papers; one BioProject split into two, since it
bundles two distinct procedures under one accession):

| Cohort | Accession | Design | Subjects | Notes |
|---|---|---|---|---|
| Cholecystectomy | PRJNA1188648 | Baseline/6M/12M | 10 | Cleanly auto-parseable from SRA `sample_title` |
| IleocecalResection | PRJNA514452 | 0/1/3/6 months | 54 | Fecal-only subset (189/509 runs) — biopsy-specimen runs excluded, not comparable to other cohorts' stool samples |
| Ileostomy | PRJNA1480144 (IL arm) | Pre/Post | 44 | Only 21/44 subjects have both timepoints |
| SDT | PRJNA1480144 (SD arm) | Pre/Post | 34 | Stent-based diversion technique; only 11/34 subjects have both timepoints |

### How to run

Same 8-stage pipeline as phase 1, extended with array tasks 5-8 (the 4 new
cohorts) in `02-05_*.sbatch`. `06_analysis.sbatch` takes a `PHASE2_COHORTS`
flag — default `0` behaves exactly as the phase-1-only pipeline always has;
set to `1` to include all 8 cohorts through the generalized
`combineCountTables.R`/`compareStudies_16S.R`/`Heatmap.R`/`compareStudies_SV.R`:

```
PHASE2_COHORTS=1 sbatch scripts/06_analysis.sbatch
```

## Workflow (phase 3 — bariatric positive control)

Phase 2 tested whether the signature generalizes to *other* GI surgeries.
Phase 3 asks a narrower, more targeted question: does it generalize to
*more bariatric surgery*, and specifically, is it tied to RYGB or to
bariatric surgery generally? The original 4 phase-1 cohorts are not
generically "bariatric" — the authors' own paper states all 4 underwent
**RYGB specifically**. "Bariatric surgery" is an umbrella term covering
RYGB, sleeve gastrectomy, gastric banding, and duodenal switch, so a
candidate matching the umbrella term isn't automatically matching the
actual procedure. Phase 3 is split into two sub-phases that ask
deliberately separated questions:

- **Phase 3b — RYGB-specific (the true positive control):** does another
  independent *RYGB* cohort share the signature?
- **Phase 3a — general bariatric, non-RYGB:** does the signature extend to
  *other* bariatric procedures (here: sleeve gastrectomy)?

This reuses phase 2's infrastructure almost entirely —
`09_find_datasets.sbatch`'s scan, `10_inspect_candidates.sh`'s readiness
check, the `cohorts/<Study>/metaData.txt` pattern, the per-cohort R script
templates, and the generalized shared scripts.

One BioProject covers both sub-phases at once, same pattern as phase 2's
`PRJNA1480144`→Ileostomy/SDT split: `PRJNA859272` (Lazaro et al. 2025,
*Obesity Surgery*) is a head-to-head SG-vs-RYGB comparative study, split
into two cohorts at integration time.

### Cohorts used

| Cohort | Sub-phase | Accession | Procedure | Subjects | Runs |
|---|---|---|---|---|---|
| RYGB_SaoPaulo | 3b | PRJNA1364303 | RYGB ± exercise | 31 | 87 |
| Lazaro_GB | 3b | PRJNA859272 (GB arm) | RYGB | 10 | 20 |
| Lazaro_SG | 3a | PRJNA859272 (SG arm) | Sleeve gastrectomy | 14 | 28 |
| Chaudhari_SG | 3a | PRJNA635168 | Sleeve gastrectomy | 17 | 34 |

Both sub-phases targeted 2 cohorts each (revised down from an original
target of 3 after an exhaustive search found no additional candidate that
didn't require either a manual/non-automatable retrieval step or came with
caveats serious enough to weaken rather than strengthen the evidence — see
[results/README.md](results/README.md)'s dataset search log for the full
candidate list and why each one was set aside).

### How to run

`scripts/02-05_*.sbatch` extended with array tasks 9-12 for the 4 phase-3
cohorts. `06_analysis.sbatch` takes a `PHASE3_COHORTS` selector (`3a` / `3b`
/ `all`, default `0`), combined with the existing `PHASE2_STUDIES` env var
(read directly by `combineCountTables.R`/`compareStudies_16S.R`/`Heatmap.R`)
to control exactly which cohorts feed the downstream comparison:

```
# OG + 3a (non-RYGB bariatric)
PHASE2_COHORTS=0 PHASE3_COHORTS=3a \
  PHASE2_STUDIES="BS,Assal,Ilhan,Afshar,Chaudhari_SG,Lazaro_SG" \
  sbatch scripts/06_analysis.sbatch

# OG + 3b (RYGB-specific)
PHASE2_COHORTS=0 PHASE3_COHORTS=3b \
  PHASE2_STUDIES="BS,Assal,Ilhan,Afshar,RYGB_SaoPaulo,Lazaro_GB" \
  sbatch scripts/06_analysis.sbatch

# OG + all 4 phase-3 cohorts together
PHASE2_COHORTS=0 PHASE3_COHORTS=all \
  PHASE2_STUDIES="BS,Assal,Ilhan,Afshar,RYGB_SaoPaulo,Lazaro_GB,Lazaro_SG,Chaudhari_SG" \
  sbatch scripts/06_analysis.sbatch

# All 12 cohorts (OG + phase 2 + phase 3) in one run
PHASE2_COHORTS=1 PHASE3_COHORTS=all sbatch scripts/06_analysis.sbatch
```

The OG-only baseline needs no new run — `PHASE2_COHORTS=0 PHASE3_COHORTS=0`
(the default) reproduces exactly the phase-1 pipeline.

**Note on parallel runs:** each `06_analysis.sbatch` invocation scaffolds
into a shared scratch pipeline directory (`$DATA_ROOT/pipeline`) and wipes
it at the start of the run. Submitting two of these runs back-to-back is
safe sequentially, but a second run's scaffold step will delete the first
run's raw output before it can be archived — either archive a completed
run's output before submitting the next one, or use a SLURM
`--dependency=afterok:<jobid>` chain and check completion before the
dependent job starts wiping the directory.

## Reporting — pairwise comparison heatmaps and scatterplots

Three small R scripts in `analysis/RScripts/` build a combined
heatmap+scatterplot report across any set of completed runs, without
re-running any statistics:

- **`pairwiseHeatmap.R <Genus_PairwiseComparison.txt> <out_prefix>`** — reads
  a run's pairwise comparison table and renders a study-timepoint ×
  study-timepoint heatmap (cell = Spearman r + BH-adjusted p), plus r/p CSV
  matrices. Same-study pairs are shaded separately since they're excluded
  from the headline correlation metric.
- **`crossGroupScatter.R <CompareStudies16S module dir> <groupA> <groupB> <title> <out_prefix>`**
  — re-plots only the scatterplots between two named groups of studies
  (e.g. "new cohorts" vs. "OG"), reusing the already-gathered
  `MixedLinearModelResults` files and `compareStudies()`/
  `plotPairwiseStudies()` from `functions.R` — a filtered re-plot, not a
  new statistical computation.
- **`titlePage.R <out.pdf> <title> <subtitle lines...>`** — a one-page
  section divider for assembling multiple PDFs into one report (merge with
  `gs -sDEVICE=pdfwrite -dBATCH -sOutputFile=out.pdf page1.pdf page2.pdf ...`).

`results/phase3/report/pairwise_comparison_report.pdf` is the combined
report built this way, covering every phase plus a taxon-level "top movers"
analysis — see [results/README.md](results/README.md) for what it found.
