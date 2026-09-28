# Gut Microbiome Signature Correlation Analysis

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
├── results/phase2/            # Extension outputs, all 8 cohorts (see "Results (phase 2)" below)
│   ├── figures/                  # Figures from the 8-cohort run (13 PDFs)
│   ├── pipeline/                 # Per-stage outputs, numbered 00-40 (gaps = skipped
│   │                              #   metagenomics modules) -- .RData LASSO model objects
│   │                              #   excluded (~635MB, rebuildable, not source work)
│   ├── predictions_lasso.tsv     # Cross-study LASSO predictions, all 8 studies
│   └── *.out                     # SLURM job logs for every successful stage run
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

## Phase 3 (proposed) — a positive control: more bariatric cohorts

Phase 2 tested whether the bariatric-surgery microbial signature
generalizes to *other* GI surgeries and found that it largely doesn't
(r drops from ~0.40 to 0.09 to -0.08 per new cohort — see "Results (phase
2)" above). That result has an alternative, less interesting explanation
worth ruling out directly: maybe cross-study correlation just weakens
whenever *any* additional, independently-collected cohort gets pooled in,
regardless of whether the surgery type matches — noisier real-world data,
different labs/protocols, more heterogeneity in general.

**The proposed test:** find additional 16S cohorts that match the
original 4 as closely as possible — a genuine pre-surgical baseline plus
≥1 post-surgical timepoint, repeat-sampled per subject — and run them
through the exact same pipeline. If these cohorts correlate with the
existing 4 at something close to the original ~0.40 (unlike phase 2's
non-bariatric cohorts), that's direct, positive evidence the signature is
real, not an artifact of "more pooled data always looks worse."

**Procedure check (2026-09-28) matters here.** The original 4 phase-1
cohorts are not generically "bariatric" — the authors' own paper states
all 4 underwent **RYGB specifically** ("three smaller publicly
available... datasets from patients who also underwent RYGB surgery"),
confirmed per cohort against each source paper (BS, Assal, Ilhan, Afshar
all RYGB). "Bariatric surgery" is an umbrella term — RYGB is one specific
procedure under it, alongside sleeve gastrectomy, gastric banding, and
duodenal switch — so a candidate matching the *umbrella* term isn't
automatically matching the *actual* procedure. This is why phase 3 is now
split into two sub-phases that ask two different, deliberately separated
questions:

- **Phase 3b — RYGB-specific (the true positive control):** does another
  independent *RYGB* cohort share the signature? A weak result here would
  be hard to explain away — this is the cleanest possible test of phase 1's
  actual finding.
- **Phase 3a — general bariatric, non-RYGB:** does the signature extend to
  *other* bariatric procedures (sleeve gastrectomy, gastric banding,
  duodenal switch, etc.)? A weak result here is ambiguous — it could mean
  "the signature isn't RYGB-specific" (still supports bariatric-specificity
  generally) or "the signature is genuinely RYGB-specific and these
  procedures differ" (a real procedure-type effect, not disconfirmation of
  phase 1). Don't treat 3a and 3b results as equivalent evidence.

Each sub-phase targets 2 papers (revised down from an original target of 3,
after two full search passes — see the "considered, not used" notes below
each table for why the natural 3rd candidates were set aside rather than
force-included).

**This reuses phase 2's infrastructure almost entirely** —
`09_find_datasets.sbatch`'s scan, `10_inspect_candidates.sh`'s
readiness check, the `cohorts/<Study>/metaData.txt` pattern, the
per-cohort R script templates (Ilhan/Afshar-style depending on timepoint
shape), and the now-generalized shared scripts
(`combineCountTables.R`/`compareStudies_16S.R`/`Heatmap.R`/
`compareStudies_SV.R`, already proven to scale past 4 studies). Phase 3 is
primarily a dataset-sourcing and onboarding task, not new engineering.

**One BioProject covers both sub-phases at once, same pattern as phase 2's
`PRJNA1480144`→Ileostomy/SDT split.** `PRJNA859272` (Lazaro et al. 2025,
*Obesity Surgery* — "Sleeve Gastrectomy and Gastric Bypass Impact in
Patient's Metabolic, Gut Microbiome, and Immuno-inflammatory Profiles") is
a head-to-head SG-vs-RYGB comparative study, 24 patients total, each
sampled at 0M (1wk pre-op) and 6M post-op. Confirmed cleanly resolvable at
the sample level — `sample_title` directly encodes subject+timepoint+
procedure (e.g. `01BP`=post-RYGB, `04S`=post-SG) — 14 SG + 10 RYGB post-op
samples, matching the paper's reported split exactly when independently
counted from the ENA filereport, not just trusted from the paper's
aggregate n. Split into two cohorts at integration time, same as phase 2.

### Phase 3b — RYGB-specific: **target met, 2/2**

| Accession | Paper | Paired subjects | Read length |
|---|---|---|---|
| `PRJNA1364303` ✅ | São Paulo, RYGB + exercise (32 women, 16 RYGB-only + 16 RYGB+exercise) | **30/31** (23 with all 3 timepoints); mixed with a mouse arm in the same BioProject — must filter on `host=Homo sapiens` AND alias containing `RYGB`, not run count alone | 285bp paired-end MiSeq → truncLen 200 |
| `PRJNA859272` (GB arm) ✅ | Lazaro et al. 2025, *Obesity Surgery* | **10/10 paired** (0M + 6M), procedure confirmed per-sample via `sample_title` suffix "BP" | ~245bp ×2 MiSeq V4 → truncLen 200 |

**Considered, not used:** `PRJNA1227268` ("Examining spatial microbiome
variations across GI tract regions in obesity," *Sci Rep* 2025) — real,
public, confirmed RYGB via explicit paper text, but set aside rather than
included as a 3rd entry: only 9/51 subjects paired, and its post-op sample
was collected **days 4-5 after surgery** — every other cohort in this
entire project measures weeks-to-months out (BS is the fastest at 1
month), so whether a surgical signature is even detectable that early is
a genuinely open question. Adding it would have diluted rather than
strengthened this sub-phase's evidence.

### Phase 3a — general bariatric, non-RYGB: **target met, 2/2**

| Accession | Paper | Procedure | Paired subjects | Read length |
|---|---|---|---|---|
| `PRJNA635168` ✅ | Chaudhari et al. 2021, *Cell Host & Microbe* — gut-liver axis | Sleeve gastrectomy | **17/17 (100%)** — cleanest pairing found; `sample_title` is literally `<subjectID>-Pre-SG`/`-Post-SG` | 250bp paired-end MiSeq → truncLen 200 |
| `PRJNA859272` (SG arm) ✅ | Lazaro et al. 2025, *Obesity Surgery* (same paper as the 3b entry above) | Sleeve gastrectomy | **14/14 paired** (0M + 6M), procedure confirmed per-sample via `sample_title` suffix "S" | ~245bp ×2 MiSeq V4 → truncLen 200 |

**Considered, not used:** `PRJNA1133160` (One-Anastomosis Gastric Bypass,
La Trobe University pilot) — technically clean (7/7 paired, cleanly
parseable metadata, public SRA, comfortably long reads) but set aside for
two stacked reasons rather than one: no peer-reviewed paper exists for
this specific dataset (only the depositing lab's own BioProject
description and an unpublished thesis draft back the procedure/design
claims), and n=7 is far smaller than every other cohort in this project.
Either alone would be a documentable caveat; together they made this a
meaningfully weaker addition rather than a genuine 3rd data point. Worth
revisiting if the underlying thesis is ever published.

**`OMIX009130` fully investigated and ruled out (2026-09-28)** — not just a
tooling gap, a hard dead end for two independent reasons, confirmed
directly against the live OMIX entry:
1. It's marked **controlled access** under China's Ministry of Science and
   Technology Human Genetic Resources framework — a formal international-
   collaboration approval process, not a simple API key or FTP credential.
2. Even past that gate, **no raw FASTQ files are deposited under this
   accession at all** — only an already-processed ASV taxonomy table
   (`.xls`, 8.8MB), which wouldn't feed this repo's raw-reads pipeline
   (`02_download.sbatch` → `03_dada2.sbatch`) regardless of access. The
   page itself states the data "has not yet been registered" for release.

No download-path engineering would fix this — the raw reads aren't
retrievable by any automated mechanism, registered or not. Not
recommended for further pursuit unless the paper's authors are contacted
directly for raw FASTQ files.

Gastric banding and duodenal switch single-procedure public cohorts were
searched for across three separate passes and consistently came up empty —
every human duodenal switch/BPD 16S study found was either not public,
cross-sectional, single-timepoint (no baseline), hosted on a non-SRA/ENA
repository, or blocked by the same unresolvable-procedure-mix problem
below. This appears to be a genuine gap in what's publicly available, not
a search-effort gap. (One omega-loop/OAGB candidate, `PRJNA1133160`, was
found — see "considered, not used" above.)

**Bottom line: both sub-phases are complete at 2/2 each.** The target was
revised down from 3 to 2 per sub-phase after two full search passes found
no additional candidate that didn't require either a manual/non-automatable
retrieval step or came with stacked caveats serious enough to weaken rather
than strengthen the evidence (see "considered, not used" above). Not
pursuing `PRJEB28869` further (below) — its supplement requires a manual
human-browser download, which wasn't pursued.

### Doesn't fit either bucket

- `PRJEB28869` (Shen et al. 2019, *SOARD*, DOI 10.1016/j.soard.2019.05.038)
  — has an explicit *Supplemental Table 3: "Comparison between RYGB and SG
  surgical groups in the Barcelona cohort"* that would almost certainly
  resolve its procedure-mix problem, but automated retrieval failed (an
  NIH-manuscript deposit gated behind a JS proof-of-work challenge that
  `curl`/`WebFetch` can't solve; ScienceDirect/ResearchGate both return
  HTTP 403). Not pursued further — would require a human manually
  downloading `NIHMS1531476-supplement-1.pdf` from
  `https://pmc.ncbi.nlm.nih.gov/articles/PMC6722012/`, and even then, ENA's
  own metadata for this BioProject is separately messy (bare sequential
  `sample_alias`, incorrect `geo_loc_name` on every sample) and would need
  manual cross-referencing by run order rather than a clean field match.
- `PRJEB39382` (Ben Izhak et al. 2021, *mSystems*) — previously listed as a
  top candidate (30/86 paired, largest paired-n found at the time), but its
  patients are an unresolvable ~1/3-each mix of Sleeve/Omega Loop/RYGB with
  no per-subject procedure labels available. Set aside unless the paper's
  authors can be reached for a per-subject procedure breakdown.
- **`PRJNA639545`** (*Gut Microbes* 2022) — **ruled out, not just blocked.**
  Structurally the best timepoint design found in this whole search (26
  RYGB + 14 SG, 4 clean timepoints: baseline/post-diet/1mo/3mo), but
  confirmed via the actual supplementary files (this journal is genuinely
  open access, retrieved via Europe PMC's API) that the supplement contains
  only 3 figures, no tables. The main text's Table 1 gives per-patient
  Sex/Age (e.g. "F/58") with no numeric ID, and BioSample XML carries no
  sex/age fields to cross-match against — there is no crosswalk anywhere
  accessible. A genuine dead end, not just an unretrieved one.
- `PRJEB48942` (Han et al. 2022, *Diabetes, Obesity and Metabolism*) — the
  richest timepoint ladder of anything found (preop/1mo/3mo/6mo/12mo, 37/41
  paired, 21 with all 5 timepoints), procedure still not confirmed against
  the paper. Read length also borderline (~150bp mean, same tier as
  Ilhan's own truncLen) — needs a real per-FASTQ check either way.

**Ruled out entirely** (mixed human+mouse translational studies,
specimen-purity failures, cross-sectional designs, or non-16S data —
corrects earlier versions of this list):
- `PRJNA727576` — only 6 of 358 "runs" are actually human; the rest are
  gnotobiotic mice inheriting the human donor's metadata tags. Too small.
- `PRJNA941111` — `isolation_source` is 100% gastric tissue (biopsy), not
  fecal. Same specimen-purity failure as IleocecalResection's excluded
  biopsy runs in phase 2.
- `PRJNA951705` — only 8/47 runs are actually human; resolves to 4 real
  subjects, too small.
- `PRJEB88699` — confirmed ~74bp mean read length, well under even
  Ilhan's 150bp floor.
- `PRJNA1118203`, `PRJNA1199984` — no parseable subject/timepoint
  attributes anywhere in BioSample XML; would need each paper's own
  supplementary table to use at all.
- `PRJNA1093424` ("Pediatric Bariatric Surgery") — confirmed via its paper
  (*Gut Microbes* 2025) that only the murine FMT-recipient samples are 16S;
  the human samples are shotgun metagenomic.
- `PRJNA321731` — confirmed cross-sectional per the source paper (Ilhan et
  al. 2017, *ISME J*): RYGB-alumni, band-alumni, normal-weight, and
  pre-surgery-obese are four different sets of people, not the same
  subjects before/after.
- `PRJNA885527`, `PRJNA355361` — confirmed rat/mouse via ENA project XML
  `SCIENTIFIC_NAME`, not human.
- `PRJNA1233226` — confirmed to be an FMT-into-mice tumor study using
  post-surgery human stool only as a transplant source, same translational
  pattern as `PRJNA727576`.

**Methodology note for future scans:** `09_find_datasets.sbatch`/
`10_inspect_candidates.sh` currently score off `sample_title`/`sample_alias`/
`library_name`/`scientific_name` only, which several ruled-out
translational studies above sailed through (`scientific_name` says
`"metagenome"` for both human and mouse arms). A `host` BioSample-attribute
check would catch this class of false positive automatically, and neither
script currently checks the *specific surgical procedure*, nor whether an
accession is even hosted on SRA/ENA vs. another repository (OMIX/NGDC,
etc.) — worth adding all three if phase 3 scans more candidates.

Verify each with `scripts/10_inspect_candidates.sh <ACCESSION>` before
committing to one, exactly as was done for the phase-2 cohorts.
