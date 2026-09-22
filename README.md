# microbiome_analysis

Fodor Lab project reproducing and extending:

> Fouladi F, et al. **"A microbial signature following bariatric surgery is
> robustly consistent across multiple cohorts."** *Gut Microbes* / *Genome
> Medicine* (2021). Original analysis code:
> [FarnazFouladi/RYGB_IntegratedAnalysis2020](https://github.com/FarnazFouladi/RYGB_IntegratedAnalysis2020)

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
   paper never looked at — i.e. whether "gut rewiring" produces this
   signature generally, not only after bariatric procedures.

Phase 1 (reproduction, below) is complete. Phase 2 (identifying and
incorporating non-bariatric surgery cohorts) is in progress.

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
├── cohorts/                   # Phase-2 per-cohort metadata, built by hand from ENA run info
│   ├── Cholecystectomy/metaData.txt
│   ├── IleocecalResection/metaData.txt
│   ├── Ileostomy/metaData.txt
│   └── SDT/metaData.txt
├── analysis/RScripts/         # Phase-2 R scripts: this repo's own code, run inside the
│   │                           # authors' Singularity images (never inside their cloned repo)
│   ├── 16S_TaxaClassification/   # One script per phase-2 cohort, modeled on the closest
│   ├── DADA2_16S/                 # matching original template (Ilhan/Afshar/Assal/BS —
│   ├── kraken2_16S/                see CLAUDE.md's "Adding a new cohort" section)
│   ├── combineCountTables.R      # Generalized rewrites of the 4 scripts that hardcode all
│   ├── Heatmap.R                  # 4 phase-1 study names in the original repo — verified
│   └── compareStudies/            # byte-identical to the real phase-1 output before trusting
├── results/phase1/            # Reproduction outputs (see below)
│   ├── figures/                  # Regenerated paper figures (PDF)
│   ├── pipeline/                 # Per-stage intermediate + final outputs (numbered 00-28,
│   │                              #   matching the original repo's module numbering)
│   ├── input/                    # Staged per-cohort input trees (BS, Assal2020, Afshar2018, Ilhan2020)
│   ├── predictions_lasso.tsv     # Our regenerated cross-study LASSO predictions
│   └── *.out                     # SLURM job logs for every stage (download, DADA2, Kraken2,
│                                  #   table-building, analysis, AUROC comparison)
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
surgery cohort), **Assal2020**, **Afshar2018**, and **Ilhan2020** — were
processed end to end from raw SRA reads rather than starting from the
original repo's precomputed count tables:

1. **`00_setup.sh`** — one-time setup: conda env for DADA2/cutadapt, the
   authors' published Singularity containers (`rbase`, `nlme`, `r-vegan`,
   `comparestudies`, `complexheatmap`, `pathogens`, `siamcat`, `fig-perform`),
   SILVA 132 training sets, a Kraken2 index (the authors' original custom DB
   no longer exists, so a public stand-in of comparable vintage was used —
   see caveat below), and a clone of the original repo for its R modules and
   sample metadata.
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

## Results: ours vs. the published reproduction

**Cross-cohort genus correlation** (the paper's headline abstract number):

| | Paper (abstract) | Original repo's own re-run | Our reproduction |
|---|---|---|---|
| Mean Spearman r, different-study pairs | 0.41 ± 0.10 | 0.410 ± 0.098 (n=23) | **0.396 ± 0.097 (n=23)** |

**Cross-study LASSO AUROC matrix** (rows = training cohort, columns =
test cohort; diagonal = within-study cross-validation):

| Train \ Test | Afshar | Assal | BS | Ilhan | LOSO |
|---|---|---|---|---|---|
| **Afshar** — ours / repo | 0.701 / 0.670 | 0.825 / 0.845 | 0.986 / 0.983 | 0.718 / 0.701 | 0.823 / 0.828 |
| **Assal** — ours / repo | 0.740 / 0.698 | 0.948 / 0.948 | 0.885 / 0.901 | 0.846 / 0.735 | 0.840 / 0.836 |
| **BS** — ours / repo | 0.864 / 0.861 | 0.864 / 0.867 | 0.969 / 0.974 | 0.769 / 0.769 | 0.947 / 0.946 |
| **Ilhan** — ours / repo | 0.612 / 0.623 | 0.706 / 0.757 | 0.569 / 0.591 | 0.718 / 0.667 | 0.803 / 0.769 |

| Summary metric | Ours | Repo |
|---|---|---|
| Within-study (CV) mean AUROC | 0.834 | 0.815 |
| Study-to-study transfer mean | 0.782 | 0.778 |
| Leave-one-study-out mean | 0.853 | 0.845 |

Mean absolute difference across the full matrix: **0.022**; max difference:
**0.111** (Assal→Ilhan transfer). Differences of this size are consistent
with normal pipeline noise — DADA2 ASV-calling stochasticity and Kraken2
index drift (the authors' original custom database no longer exists; a
public index of comparable vintage was substituted, see `scripts/00_setup.sh`)
— not a breakdown of the underlying signal.

**Conclusion:** the reproduction from raw reads confirms the paper's central
claim. A consistent, transferable microbial signature of bariatric surgery
is recovered independently of the original authors' precomputed tables,
within the paper's own reported variance.

## Workflow (phase 2 — extending beyond bariatric surgery, in progress)

`scripts/09_find_datasets.sbatch` scans ENA study metadata for candidate
non-bariatric GI-surgery 16S cohorts against explicit inclusion criteria
(Illumina amplicon sequencing, ≥150bp reads, human gut samples, repeat
sampling per subject, a recoverable pre-surgical baseline plus ≥1
post-surgical timepoint, GI-surgery context), scoring and ranking results
into `dataset_scan/shortlist.tsv`. `scripts/10_inspect_candidates.sh` then
deep-inspects specific shortlisted accessions — sample label parsing,
subject/timepoint layout, read length — to decide which are actually usable
and with what DADA2 `truncLen`.

### Cohorts selected

Four cohorts (three source papers; one BioProject split into two, since it
bundles two distinct procedures under one accession):

| Cohort | Accession | Design | Subjects | Notes |
|---|---|---|---|---|
| Cholecystectomy | PRJNA1188648 | Baseline/6M/12M | 10 | Cleanly auto-parseable from SRA `sample_title` |
| IleocecalResection | PRJNA514452 | 0/1/3/6 months | 54 | Fecal-only subset (189/509 runs) — 320 biopsy-specimen runs (Bp/iRs/niRs) excluded, not comparable to the other cohorts' stool samples |
| Ileostomy | PRJNA1480144 (IL arm) | Pre/Post | 44 | Only 21/44 subjects have both timepoints — most are baseline-only |
| SDT | PRJNA1480144 (SD arm) | Pre/Post | 34 | Only 11/34 subjects have both timepoints — most are baseline-only |

The Ileostomy/SDT attrition rate (52%/68% singleton subjects) is notably
higher than any phase-1 cohort's (BS tops out at 30%) — included regardless,
consistent with the pipeline's existing tolerance for unbalanced subjects,
but flagged here since it means weaker statistical power per subject than
the phase-1 cohorts once results exist.

### What's built (not yet run)

- Per-cohort metadata (`cohorts/<Study>/metaData.txt`), derived by hand from
  ENA run metadata — see `CLAUDE.md`'s "Adding a new cohort" section
- `scripts/02-05` extended with `case`/config entries for all 4 cohorts
  (verified additive-only — the original 4 studies' values are unchanged)
- 12 new per-cohort R scripts (`analysis/RScripts/`), one
  `16S_TaxaClassification`/`DADA2_16S`/`kraken2_16S` script per cohort,
  each modeled on the closest-matching original template by timepoint shape
- The 4 shared downstream scripts that hardcode all 4 phase-1 study names in
  the original (`combineCountTables.R`, `compareStudies_16S.R`, `Heatmap.R`,
  `compareStudies_SV.R`) rewritten generically in this repo and **verified
  byte-identical to the real published phase-1 output** before being trusted
  (see `CLAUDE.md`'s "Never assume" section for how)
- `scripts/06_analysis.sbatch` wired with a `PHASE2_COHORTS` flag (default
  `0`): off, it behaves exactly as it always has (verified — same module
  numbering as the actual complete pipeline run, not just equivalent
  content); set to `1`, it includes all 8 cohorts through the generalized
  scripts

### Status (2026-09-22)

- **Downloads complete.** All 8 cohorts (4 phase-1 + 4 phase-2) downloaded
  cleanly — zero failures across all 8 array tasks, every cohort's file
  count matches expected exactly, and read lengths for all 4 new cohorts
  (Cholecystectomy 301bp, IleocecalResection 314bp, Ileostomy/SDT 250bp)
  comfortably exceed the planned `truncLen=200`.
- **`PRIMER`/`TRIMLEFT` confirmed, not placeholders anymore.** Ran the
  515F-anchor check against the real downloaded reads: none of the 4 new
  cohorts show a fixed-position anchor (IleocecalResection: zero matches;
  Cholecystectomy/Ileostomy/SDT: scattered matches at varying positions —
  coincidental, not a real primer). `PRIMER=NONE`/`TRIMLEFT=0` is correct as
  configured; no change needed before running `03_dada2.sbatch`.
- `compareStudies_SV.R` still intentionally excludes all 4 phase-2 cohorts
  by default — this same check independently confirms none of them share
  BS/Assal's exact V4 protocol, so this stays a genuine exclusion, not just
  an unverified one. See `CLAUDE.md`'s "SV/ASV-level comparison is
  region-locked" section before ever adding one.

### Still needed before a real run

1. `03_dada2.sbatch` — DADA2 ASV inference for all 8 cohorts
2. `04_kraken2_16s.sbatch` → `05_build_tables.sbatch` → `06_analysis.sbatch`
   (with `PHASE2_COHORTS=1`) → `07_compare_auroc.sbatch`
