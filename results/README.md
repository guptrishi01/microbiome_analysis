# Results

This file holds every numeric finding and interpretation produced by the
pipeline described in the top-level [README.md](../README.md). That file
covers what the pipeline *is* and how to run it; this one covers what running
it actually found. Each phase's raw output (figures, pipeline module
directories, SLURM logs, LASSO predictions) lives in `results/phaseN/`,
referenced below.

## Phase 1 — reproduction vs. the published paper

Pipeline output: `results/phase1/`.

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
public index of comparable vintage was substituted, see `CLAUDE.md`'s
"Cluster environment" section) — not a breakdown of the underlying signal.

**Conclusion:** the reproduction from raw reads confirms the paper's central
claim. A consistent, transferable microbial signature of bariatric surgery
is recovered independently of the original authors' precomputed tables,
within the paper's own reported variance.

### Verification notes (phase 1)

Every stage's output was checked directly against the original repo's own
published tables before being trusted, not assumed equivalent from reading
the code:

- **Kraken2 database substitution** is the single biggest source of the
  small numeric differences above — the authors' original custom database
  no longer exists, so `00_setup.sh` substitutes a public index of
  comparable vintage. Expected and documented, not chased to zero.
- Read-length/primer settings (`truncLen` 200bp for BS/Assal/Afshar, 150bp
  for Ilhan) reproduce the original repo's own ASV length distributions
  exactly.

## Phase 2 — extension to non-bariatric GI surgery

Pipeline output: `results/phase2/`. Cohorts used: Cholecystectomy,
IleocecalResection, Ileostomy, SDT — see the top-level README's Phase 2
section for accessions/design/subject counts.

### Pipeline run notes

- **Downloads:** all 8 cohorts (4 original + 4 new), zero failures, every
  file count matched expected exactly. Read lengths for the 4 new cohorts
  (Cholecystectomy 301bp, IleocecalResection 314bp, Ileostomy/SDT 250bp)
  comfortably exceeded `truncLen=200`.
- **Primer/trim settings confirmed against real reads**, not left as
  placeholders: none of the 4 new cohorts showed a fixed-position 515F
  anchor, confirming `PRIMER=NONE`/`TRIMLEFT=0` was correct. This also
  confirms none of them share BS/Assal's exact V4 protocol, so
  `compareStudies_SV.R` correctly excludes all 4 by default (see
  `CLAUDE.md`'s "SV/ASV-level comparison is region-locked").
- **DADA2:** all 8 succeeded, every ASV length matched its configured
  `EXPECT_LEN`. Two samples dropped below 1000 reads in IleocecalResection
  (`SRR13357992`, `SRR13358126`) — caught by the script's own printed
  warning; excluded later by `norm()`'s 1000-read filter (~1% of that
  cohort's 189 samples).
- **Lower chimera-removal read retention** for Cholecystectomy (71.5%),
  Ileostomy (67.1%), SDT (66.9%) vs. 90-99% for every other cohort. Every
  sample in all 3 still had tens of thousands of reads of margin above the
  1000-read cutoff, so no sample-dropout risk — but this does mean reduced
  per-taxon estimation precision for these 3 cohorts specifically, a real
  interpretation caveat for their weaker correlation numbers below.
- **Kraken2:** all 8 succeeded, every mpa-report count matched expected.
  Three cohorts (Afshar, Ileostomy, SDT) each had one isolated
  low-classification-rate sample, confirmed as single outliers.
- **Table-building:** all 8 succeeded; genus-level correlation against the
  original repo's published tables for the 4 phase-1 cohorts was r=0.82-0.85,
  consistent with the Kraken2-database-substitution pattern, not a new
  discrepancy. One real bug caught and fixed here: the first submission
  crashed on a missing `COHORTS` argument to `build_tables.py` — fixed,
  verified the regenerated phase-1 tables were still byte-identical to the
  already-published output, then re-ran.
- **`06_analysis.sbatch`:** first submission crashed — IleocecalResection's
  TaxaClassification module hit a row-count mismatch (`norm()` drops 2
  sub-1000-read samples, but the original `cbind(t1_norm, meta)` pattern
  assumes row counts always match). Fixed in all 4 new cohort scripts
  (subset `meta` to `rownames(t1_norm)` immediately before each `cbind`),
  verified against the real crashed data (187/187 rows, consistent across
  all 6 ranks) and against a cohort where nothing gets dropped (byte-identical
  output before/after the fix, confirming a true no-op elsewhere).
  Resubmitted, completed clean: all 34 modules succeeded.

### Genus-level cross-study correlation

The full 8-cohort, 105-pair comparison came back at **r = 0.191 ± 0.224**,
far below the paper's number. Before treating that as a finding, this was
checked for a pipeline bug: the 23 phase-1-only pairs *within* this same
105-pair table gave **r = 0.397 ± 0.097** — essentially identical to the
already-verified phase-1-only result — confirming the generalized
`compareStudies_16S.R` is correct and the drop comes entirely from the new
cohorts, not from anything broken in the shared code.

Breaking the new cohorts' cross-study correlation down individually:

| Cohort | r (different-study pairs) |
|---|---|
| IleocecalResection | 0.188 ± 0.113 (n=30) — weak positive |
| Ileostomy | 0.093 ± 0.136 (n=9) — near zero |
| Cholecystectomy | -0.010 ± 0.097 (n=26) — no signal |
| SDT | -0.078 ± 0.065 (n=8) — slightly negative |

### Cross-study LASSO AUROC

Within-study performance varies sharply by cohort — Afshar 0.717, Assal
0.954, BS 0.984, Cholecystectomy 0.745, **IleocecalResection 0.533** (barely
above chance), Ileostomy 0.892, Ilhan 0.744, SDT 0.783 — and pooled
cross-study transfer drops from the phase-1-only run's 0.778 to **0.619**,
with LOSO mean dropping from 0.845 to **0.683**.

Pooling the noisier phase-2 cohorts into the same LASSO training set also
measurably hurt the phase-1 cohorts' own LOSO performance, not just the new
cohorts' transfer numbers:

| Phase-1 cohort | LOSO AUROC, phase-1-only | LOSO AUROC, pooled with all 8 | Δ |
|---|---|---|---|
| Afshar | 0.828 | 0.659 | -0.169 |
| Assal | 0.836 | 0.813 | -0.023 |
| BS | 0.946 | 0.585 | -0.361 |
| Ilhan | 0.769 | 0.650 | -0.120 |

### Interpretation (phase 2)

The bariatric-surgery microbial signature does not straightforwardly
generalize to the other GI-surgery types tested here — IleocecalResection
shows a real but much weaker signal, and Cholecystectomy/SDT show
essentially none. IleocecalResection's own within-study AUROC (0.533) being
little better than chance, despite comparable data quality to the phase-1
cohorts (90.6% chimera retention, a complete repeat-sampling design), is
evidence this reflects a real procedural/biological difference rather than
only the data-quality caveats above — a fundamentally different patient
population (Crohn's disease) and a 4-level, gradual timepoint design
plausibly make the surgical signature itself harder to detect, not just
harder to measure. That said, the documented noise caveats for
Cholecystectomy/Ileostomy/SDT (lower chimera retention, mostly-unpaired
subjects) mean their near-zero numbers shouldn't be read as definitive proof
of *no* signal, only that none was detected at this sample size and data
quality.

## Phase 3 — positive control: more bariatric cohorts

Pipeline output: `results/phase3/OG_plus_3a/`, `results/phase3/OG_plus_3b/`,
`results/phase3/all12/`. Cohorts used: RYGB_SaoPaulo, Lazaro_GB (3b,
RYGB-specific), Lazaro_SG, Chaudhari_SG (3a, non-RYGB bariatric) — see the
top-level README's Phase 3 section for accessions/design/subject counts and
the full dataset-search log.

### Onboarding verification

- **Dataset sizes:** RYGB_SaoPaulo 31 subjects/87 runs, Lazaro_GB 10/20,
  Lazaro_SG 14/28, Chaudhari_SG 17/34. Phase 3b total 41 subjects, Phase 3a
  total 31 subjects.
- Two real bugs caught and fixed before any job was submitted: a doubled
  `RYGB_RYGB_SaoPaulo_` Kraken2 file prefix, and a missing `prepost`
  metadata column in RYGB_SaoPaulo's TaxaClassification script that
  `combineCountTables.R` requires — both would have caused a silent failure
  or crash downstream, caught by review before the first submission.
- **Primer check found a real, mixed result** — not the uniform "no primer"
  pattern phase 2 had. RYGB_SaoPaulo and Chaudhari_SG showed the same
  phase-2 signature (no consistent anchor) — `PRIMER=NONE` confirmed
  correct. **Lazaro_GB and Lazaro_SG did not**: ~90% of their first 5000
  reads carried the exact, untrimmed degenerate 515F primer
  (`GTGYCAGCMGCCGCGGTAA`) — the same situation as BS16S/Assal in phase 1.
  The first submission for these two cohorts ran with the placeholder
  `PRIMER=NONE` and was discarded and rerun with
  `PRIMER=GTGYCAGCMGCCGCGGTAA` (cutadapt-based removal) — confirmed correct
  by the rerun's ASV length distribution collapsing to a clean, single
  200bp peak (it was not before the fix).
- DADA2: all 4 succeeded after the primer fix, every ASV length matched
  `EXPECT_LEN=200`, no samples dropped below the 1000-read threshold.
  Kraken2: all 4 succeeded, classification rates 98.3-100%.
- Both `OG_plus_3a` and `OG_plus_3b` runs completed 28/28 modules ok.
  Pair counts in each run's `Genus_PairwiseComparison.txt` matched the
  expected C(n,2) for that sub-phase's total timepoint count (45 pairs / 10
  timepoints for 3a; 55 pairs / 11 timepoints for 3b).

### Genus-level cross-study correlation

| Configuration | Cohorts added | r (different-study pairs) | n pairs |
|---|---|---|---|
| **OG** (phase-1 baseline) | — | **0.396 ± 0.097** | 23 |
| **OG+3a** (non-RYGB bariatric) | Chaudhari_SG, Lazaro_SG (both SG) | **0.313 ± 0.158** | 40 |
| **OG+3b** (RYGB-specific) | RYGB_SaoPaulo, Lazaro_GB (both RYGB) | **0.431 ± 0.110** | 49 |
| **All 12 cohorts** (OG + phase 2 + phase 3) | all 8 added cohorts | **0.192 ± 0.206** | 180 |

### Interpretation (phase 3)

Adding two more **RYGB** cohorts (3b) does not weaken the signature at all —
it comes back essentially at, or slightly above, the phase-1 baseline
(0.431 vs. 0.396). Adding two more **bariatric-but-non-RYGB** cohorts (3a)
does measurably weaken it (0.313), a real directional decline — but nowhere
close to phase 2's collapse when genuinely unrelated GI surgeries were
pooled in (0.191). Pooling all 12 cohorts together collapses the
correlation back down to 0.192, essentially matching phase 2 alone — the 4
unrelated-surgery cohorts dilute the signal even with the RYGB-positive
cohorts folded in.

Read together, phases 2 and 3 triangulate the same conclusion from multiple
directions: the signature doesn't generalize to unrelated GI surgery; it
*does* generalize cleanly to more RYGB cohorts; and it generalizes only
partially to other bariatric procedures. That pattern — strong
within-procedure transfer, partial cross-procedure transfer, no
cross-surgery-type transfer — is most consistent with the signature being
substantially **RYGB-specific** rather than a general "any bariatric
surgery" or "any gut surgery" effect, though 3a's n=2 cohorts (both sleeve
gastrectomy specifically, not gastric banding or duodenal switch) means
this is evidence toward that conclusion, not a final word on
bariatric-procedure specificity as a whole.

## Top movers and antibiotic cross-reference

Full report (heatmaps + scatterplots for every phase, plus this analysis in
detail): `results/phase3/report/pairwise_comparison_report.pdf` (33 pages).
Underlying r/p matrices: `results/phase3/report/*.csv`.

Across all 12 cohorts' per-taxon mixed-linear-model results (BH-adjusted
p < 0.05 for "significant"):

- **Veillonella** is the single most consistent post-surgery increaser —
  significantly increased in 10 of 12 cohort-timepoints, spanning OG,
  Phase 2, *and* Phase 3b (mean effect size s=1.26). It also produces two of
  the three largest single increases found anywhere (RYGB_SaoPaulo T2/T3,
  s=1.97/1.87).
- **Streptococcus** (and several others — Rothia, Fusobacterium,
  Granulicatella, Actinomyces, Haemophilus, Klebsiella, Escherichia/Shigella)
  are the cleanest **RYGB-specific** consistent increasers: significant only
  in OG and/or Phase 3b cohorts, never in any Phase 2 or Phase 3a cohort.

**Antibiotic cross-reference.** Per-dataset antibiotic use, pulled from each
source paper:

| Cohort | Antibiotic use |
|---|---|
| BS | Prophylactic antibiotics routinely given before surgery; drug/dose not named in the paper |
| Assal, Afshar, Ilhan, RYGB_SaoPaulo, Chaudhari_SG | Not reported or not accessible |
| Lazaro_GB, Lazaro_SG | Explicitly **no** antibiotic prophylaxis used |

Veillonella still rises significantly in Lazaro_GB (s=1.63, adj.p=0.041) —
a cohort with **no antibiotic exposure at all**. The same taxon dominating
the top-movers list in antibiotic-exposed cohorts also rises in a cohort
with none, arguing against antibiotics being the primary driver — this
looks like a genuine post-surgical effect (altered gastric acid, bile flow,
oxygen gradient), not purely an antibiotic-selection artifact.

Caveat: BS is the one cohort with *confirmed* antibiotic use and shows the
strongest Klebsiella/Veillonella/Streptococcus increases of any cohort —
antibiotic contribution there specifically cannot be ruled out, since the
drug was never named. The original Fouladi et al. paper's own Discussion
raises exactly this possibility ("routinely administered operative
prophylactic antibiotics... might be the reason for the increase in
opportunistic pathogens after surgery"). No source paper across all 12
cohorts names a specific antibiotic drug, so any taxon-to-drug mechanistic
link beyond this remains unconfirmed, general pharmacology context rather
than a verified per-study match.

## Dataset search log (phase 3 candidates not used)

A record of every candidate dataset investigated for phase 3 beyond the 4
actually used, and why each was set aside — kept for anyone extending the
search later.

**Phase 3b (RYGB-specific) — considered, not used:**
`PRJNA1227268` ("Examining spatial microbiome variations across GI tract
regions in obesity," *Sci Rep* 2025) — real, public, confirmed RYGB via
explicit paper text, but set aside: only 9/51 subjects paired, and its
post-op sample was collected days 4-5 after surgery — every other cohort in
this project measures weeks-to-months out, so whether a surgical signature
is even detectable that early is a genuinely open question.

**Phase 3a (non-RYGB bariatric) — considered, not used:**
`PRJNA1133160` (One-Anastomosis Gastric Bypass, La Trobe University pilot)
— technically clean (7/7 paired) but no peer-reviewed paper exists for this
specific dataset, and n=7 is far smaller than every other cohort in this
project. Worth revisiting if the underlying thesis is ever published.

**`OMIX009130`** — fully investigated and ruled out: marked controlled
access under China's Ministry of Science and Technology Human Genetic
Resources framework, and even past that gate, no raw FASTQ files are
deposited under this accession at all (only an already-processed ASV
taxonomy table). Not retrievable by any automated mechanism.

Gastric banding and duodenal switch single-procedure public cohorts were
searched for across three separate passes and consistently came up empty —
every human duodenal switch/BPD 16S study found was either not public,
cross-sectional, single-timepoint, hosted on a non-SRA/ENA repository, or
blocked by an unresolvable procedure-mix problem. This appears to be a
genuine gap in what's publicly available, not a search-effort gap.

**Doesn't fit either bucket:**
- `PRJEB28869` (Shen et al. 2019, *SOARD*) — has a supplement that would
  likely resolve its procedure-mix problem, but automated retrieval failed
  (JS proof-of-work challenge `curl`/`WebFetch` can't solve; ScienceDirect/
  ResearchGate return HTTP 403). Not pursued — would need a human manually
  downloading the supplement.
- `PRJEB39382` (Ben Izhak et al. 2021, *mSystems*) — largest paired-n found
  (30/86), but patients are an unresolvable ~1/3-each mix of
  Sleeve/Omega Loop/RYGB with no per-subject procedure labels.
- **`PRJNA639545`** (*Gut Microbes* 2022) — ruled out, not just blocked.
  Structurally the best timepoint design found (26 RYGB + 14 SG, 4 clean
  timepoints), but confirmed via the actual (open-access) supplementary
  files that no subject-ID crosswalk exists anywhere accessible.
- `PRJEB48942` (Han et al. 2022, *Diabetes, Obesity and Metabolism*) — the
  richest timepoint ladder found (preop/1mo/3mo/6mo/12mo, 37/41 paired),
  procedure not confirmed against the paper, read length borderline (~150bp).

**Ruled out entirely** (mixed human+mouse translational studies,
specimen-purity failures, cross-sectional designs, or non-16S data):
`PRJNA727576` (6/358 runs human, rest gnotobiotic mice), `PRJNA941111`
(100% gastric-tissue biopsy, not fecal), `PRJNA951705` (8/47 runs human),
`PRJEB88699` (~74bp mean read length), `PRJNA1118203`/`PRJNA1199984` (no
parseable subject/timepoint attributes), `PRJNA1093424` (human samples are
shotgun metagenomic, not 16S), `PRJNA321731` (cross-sectional, not
repeat-sampled), `PRJNA885527`/`PRJNA355361` (rat/mouse), `PRJNA1233226`
(FMT-into-mice tumor study).

**Methodology note for future scans:** `scripts/09_find_datasets.sbatch`/
`10_inspect_candidates.sh` currently score off `sample_title`/`sample_alias`/
`library_name`/`scientific_name` only, which several ruled-out translational
studies above sailed through (`scientific_name` says "metagenome" for both
human and mouse arms). A `host` BioSample-attribute check would catch this
class of false positive automatically, and neither script currently checks
the *specific surgical procedure*, nor whether an accession is even hosted
on SRA/ENA vs. another repository (OMIX/NGDC, etc.) — worth adding all
three if the search resumes.
