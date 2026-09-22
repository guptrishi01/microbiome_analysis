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

`02_download` → `03_dada2` → `04_kraken2_16s` → `05_build_tables` →
`06_analysis` (`PHASE2_COHORTS=1`) → `07_compare_auroc` all complete for
all 8 cohorts. Phase 2's core analysis is done — see "Results" below.

- **Downloads:** all 8 cohorts, zero failures, every file count matches
  expected exactly. Read lengths for the 4 new cohorts (Cholecystectomy
  301bp, IleocecalResection 314bp, Ileostomy/SDT 250bp) comfortably exceed
  `truncLen=200`.
- **`PRIMER`/`TRIMLEFT` confirmed against real reads**, no longer
  placeholders: none of the 4 new cohorts show a fixed-position 515F
  anchor, confirming `PRIMER=NONE`/`TRIMLEFT=0` is correct. Independently
  confirms none share BS/Assal's exact V4 protocol — `compareStudies_SV.R`
  correctly excludes all 4 by default (see `CLAUDE.md`'s "SV/ASV-level
  comparison is region-locked").
- **DADA2 (`03`):** all 8 succeeded, every ASV length matches its
  configured `EXPECT_LEN`. **Two samples dropped below 1000 reads** in
  IleocecalResection (`SRR13357992`, `SRR13358126`) — the script's own
  printed warning caught these; they'll be silently excluded later by
  `norm()`'s `rowSums(table)>1000` filter (~1% of that cohort's 189
  samples, not corrected for, just noted).
- **Lower chimera-removal read retention** for Cholecystectomy (71.5%),
  Ileostomy (67.1%), SDT (66.9%) vs. 90-99% for every other cohort. Checked
  directly, not assumed: every sample in all 3 has tens of thousands of
  reads of margin above the 1000-read cutoff, so no sample-dropout risk.
  The real effect is reduced per-taxon estimation precision for these 3
  cohorts specifically — not a bug (chimera detection ran identically
  across all 8; this reflects real differences in the underlying labs'
  protocols), but a genuine interpretation caveat: if these 3 cohorts show
  weaker cross-study correlation or LASSO transfer later, that could be
  measurement noise rather than a true biological negative, and should be
  read that way rather than as a clean result.
- **Kraken2 (`04`):** all 8 succeeded, every mpa-report count matches
  expected. Three cohorts (Afshar, Ileostomy, SDT) each had one isolated
  low-classification-rate sample (80.88%/78.25%/87.79%) — confirmed single
  outliers, not systemic (the next-lowest sample in each recovers to
  91-95%+).
- **Table-building (`05`):** all 8 succeeded; genus-level correlation
  against the original repo's published tables for the 4 phase-1 cohorts
  is r=0.82-0.85, consistent with the already-documented Kraken2-database-
  substitution pattern (see "Cluster environment" in `CLAUDE.md`), not a
  new discrepancy. Hit and fixed one real bug along the way: the first
  submission crashed immediately (a missing `COHORTS` argument to
  `build_tables.py`, introduced when wiring phase-2 metadata but never
  passed at the call site) — fixed, verified the regenerated phase-1
  tables are still byte-identical to the already-published output, then
  re-ran successfully.
- **Human/Eukaryota read contamination** (present in several samples,
  expected background in 16S data) is not a residual concern — it's
  actively filtered by design: `05_build_tables.sbatch` keeps only
  `d__Bacteria`/`d__Archaea` before normalization, dropping Eukaryota/
  Viruses reads before any per-sample totals are computed. This is the one
  "noise" item that's genuinely handled by the pipeline itself, rather
  than just documented for later interpretation.

None of the above are corrupting or blocking anything — they're real,
verified properties of specific cohorts' data, not open problems, and
they're recorded here so they inform how phase-2 results get interpreted
below.

- **`06_analysis.sbatch`:** first submission crashed immediately —
  IleocecalResection's TaxaClassification module hit a row-count mismatch
  (`arguments imply differing number of rows: 187, 189`) in
  `cbind(t1_norm, meta)`. Root cause: `norm()` correctly drops the 2
  already-flagged sub-1000-read samples, but the original authors'
  `cbind(t1_norm, meta)` pattern (copied into all 4 new cohort scripts)
  assumes the row counts always match — true for all 4 phase-1 cohorts
  (none of them ever hit this), false here. Fixed in all 4 new scripts
  (subset `meta`/`meta1` to `rownames(t1_norm)` per rank, immediately
  before each `cbind`) rather than just the one that happened to trigger
  it. Verified against the real data that crashed (187/187 rows,
  consistent across all 6 ranks) and against a cohort where nothing gets
  dropped (byte-identical output before/after the fix, confirming it's a
  true no-op elsewhere). Resubmitted, completed clean: **all 34 modules
  succeeded, zero failures**, ~1h45m total runtime (`TrainModellasso`
  alone took ~63 min, training/testing LASSO across 8 studies instead of
  4).
- **`07_compare_auroc.sbatch`:** completed in 12s, no changes needed — the
  script derives its study list dynamically from `metaData_merged.txt`,
  so it handled all 8 studies without modification.

## Results (phase 2)

**Genus-level cross-study correlation** (the same metric behind the
paper's headline 0.41 ± 0.10): the full 8-cohort, 105-pair comparison came
back at **r = 0.191 ± 0.224**, far below the paper's number. Before
treating that as a finding, this was checked for a pipeline bug: the 23
phase-1-only pairs *within* this same 105-pair table give **r = 0.397 ±
0.097** — essentially identical to the already-verified phase-1-only
result (0.396 ± 0.097) — confirming `compareStudies_16S.R`'s
generalization is correct and the drop is coming entirely from the new
cohorts, not from anything broken in the shared code.

Breaking the new cohorts' cross-study correlation down individually:

| Cohort | r (different-study pairs) |
|---|---|
| IleocecalResection | 0.188 ± 0.113 (n=30) — weak positive |
| Ileostomy | 0.093 ± 0.136 (n=9) — near zero |
| Cholecystectomy | -0.010 ± 0.097 (n=26) — no signal |
| SDT | -0.078 ± 0.065 (n=8) — slightly negative |

**Cross-study LASSO AUROC** tells a consistent story, with an additional
finding of its own. Within-study performance varies sharply by cohort —
Afshar 0.717, Assal 0.954, BS 0.984, Cholecystectomy 0.745,
**IleocecalResection 0.533** (barely above chance), Ileostomy 0.892, Ilhan
0.744, SDT 0.783 — and pooled cross-study transfer drops from the
phase-1-only run's 0.778 to **0.619**, with LOSO mean dropping from 0.845
to **0.683**.

The additional finding: **pooling the noisier phase-2 cohorts into the
same LASSO training set measurably hurt the phase-1 cohorts' own LOSO
performance**, not just the new cohorts' transfer numbers:

| Phase-1 cohort | LOSO AUROC, phase-1-only | LOSO AUROC, pooled with all 8 | Δ |
|---|---|---|---|
| Afshar | 0.828 | 0.659 | -0.169 |
| Assal | 0.836 | 0.813 | -0.023 |
| BS | 0.946 | 0.585 | -0.361 |
| Ilhan | 0.769 | 0.650 | -0.120 |

**Interpretation.** The bariatric-surgery microbial signature does not
straightforwardly generalize to the other GI-surgery types tested here —
IleocecalResection shows a real but much weaker signal, and
Cholecystectomy/SDT show essentially none. IleocecalResection's own
within-study AUROC (0.533) being little better than chance, despite
comparable data quality to the phase-1 cohorts (90.6% chimera retention,
a complete repeat-sampling design), is evidence this reflects a real
procedural/biological difference rather than only the data-quality caveats
already documented above — a fundamentally different patient population
(Crohn's disease) and a 4-level, gradual timepoint design (0/1/3/6 months,
collapsed to a single pre/post split for LASSO) plausibly make the
surgical signature itself harder to detect, not just harder to measure.
That said, the documented noise caveats for Cholecystectomy/Ileostomy/SDT
(lower chimera retention, mostly-unpaired subjects) mean their near-zero
numbers shouldn't be read as definitive proof of *no* signal, only that
none was detected at this sample size and data quality.
