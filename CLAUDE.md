# CLAUDE.md

Project-specific context for this repo. See [README.md](README.md) for the paper citation, objective, results comparison, and repo structure.

## What this is

A reproduction of Fouladi et al.'s cross-cohort bariatric-surgery gut-microbiome
signature, run from raw SRA reads rather than the original repo's precomputed
tables, followed by an attempt to extend the finding to non-bariatric GI-surgery
cohorts. Phase 1 (reproduction, 4 cohorts) is done and verified against the
original repo's own published numbers. Phase 2 (extension) is in progress.

**Canonical upstream reference:**
[FarnazFouladi/RYGB_IntegratedAnalysis2020](https://github.com/FarnazFouladi/RYGB_IntegratedAnalysis2020)
— cloned read-only at `$DATA_ROOT/RYGB_IntegratedAnalysis2020`. Its original
functionality (the R analysis for the 4 phase-1 cohorts, exactly as published)
must remain intact through every phase-2 change:

- **Never write into the clone.** Confirm with `git -C
  $DATA_ROOT/RYGB_IntegratedAnalysis2020 status` before and after any session
  that touches phase 2 — it must always read "working tree clean," matching
  `origin/main`. Everything new (per-cohort metadata, new R scripts) belongs in
  *this* repo (`cohorts/`, or alongside — never inside — the clone's own tree).
- **Never change phase-1 behavior while adding phase-2 cohorts.** The 4
  original `case` branches in `02-05_*.sbatch` and the 4 original `MODULES`
  entries per stage in `06_analysis.sbatch` must stay byte-identical when
  cohorts 5+ are added — confirm with `git diff` before submitting any job, not
  just before committing. A phase-2 addition should only ever add lines, never
  edit an existing phase-1 value.

## Never assume — verify and prove it

This is the working principle for every change in this repo, not just the two
rules above. A claim that a rewrite "should behave the same" or "looks
equivalent" is not sufficient — prove it, concretely, before trusting it:

- **Prefer a direct comparison over reasoning about correctness.** When
  rewriting or extending something with known-good prior output (a script, a
  config, a data transform), run the new version against real existing data and
  diff the result against the real existing output — don't just read the code
  and conclude it's equivalent.
- **Worked example already in this repo's history:** `combineCountTables.R` was
  rewritten from 4 hardcoded study variables to a generic loop (needed to add
  phase-2 cohorts). Rather than assuming the rewrite preserved behavior, it was
  run directly (via `singularity exec .../nlme.sif Rscript`, no SLURM needed
  for something this cheap) against the real phase-1 `*TaxaClass/output/`
  directories already on disk, restricted to the original 4 studies via a
  `PHASE2_STUDIES` env var, and the output was diffed byte-for-byte against
  `results/phase1/pipeline/*_CombineCountTables*/output/`. This caught a real
  discrepancy (the `ID` column) on the first attempt — silently trusting the
  rewrite would have shipped a subtle output difference.
- **When a "test" surfaces a difference, don't just make it disappear.**
  Root-cause it first. In the case above, the difference traced to the
  *original* script's own inconsistency (it wraps `timepoint` in
  `as.character()` before combining across studies, but not `ID`, so combining
  multiple R factors via base `c()` collapsed `ID` to integer codes) — almost
  certainly an unintentional quirk in the original, but reproducing it exactly
  was the correct fix, not "improving" it, because exact reproducibility was
  the requirement.
- **No compute cluster access needed for most of this.** Singularity and every
  `.sif` image already exist locally (`images/*.sif`) — a quick correctness
  check against already-generated phase-1 output can run directly, without
  waiting on a SLURM allocation. Save `sbatch` for the real, full-scale run
  once the logic is already proven.
- **State what was actually verified, not what should logically follow from
  it.** "Untouched — confirmed via `git diff`, zero output" is a claim someone
  else can re-check. "Should be fine since nothing else changed" is not.

## Cluster environment

Two roots, set as env vars everywhere (`00_setup.sh` and every `scripts/*.sbatch`):

- `SW_ROOT` (default `$HOME/rygb`, 500G, no purge) — conda envs, Singularity images,
  reference databases (SILVA 132, Kraken2 index). Built once by `00_setup.sh`,
  idempotent (`have()` checks skip anything already present).
- `DATA_ROOT` (default `/scratch/$USER/rygb`, 5T, scratch — assume it can be purged)
  — raw reads, intermediates, pipeline outputs, and the cloned original repo
  (`$DATA_ROOT/RYGB_IntegratedAnalysis2020`). This is where the actual R modules
  and per-cohort metadata files the pipeline reads at runtime live — not in this
  git repo. `results/phase1/` in this repo is a copy of what that produced, kept
  for the record.

Modules: `mambaforge/23.11`, `singularity/4.4.1`, `sra-tools/3.1.0`,
`kraken2/2.1.3`, `python/3.12.3` (ad hoc `python3` outside a loaded module falls
back to the system default, 3.9.25 on this host — the scan scripts are stdlib-only
so both work). All compute-heavy stages are submitted via `sbatch` on the `Orion`
partition; nothing in `scripts/` should be run for real work directly on a login node.

**The authors' Kraken2 database is gone.** `00_setup.sh` substitutes a public index
of comparable vintage (`K2_URL`, default `k2_standard_20201202`) — this is the
single biggest source of the small numeric differences between our reproduction
and the original repo's published values (see README's results table). Don't
"fix" this by chasing exact number matches; the point of `07_compare_auroc.sbatch`
is confirming the signal survives the substitution, not eliminating every delta.

## Git scope

`.gitignore` excludes `conda_pkgs/`, `envs/`, `images/`, `db/`, `singularity_cache/`
(~59GB, all rebuildable via `00_setup.sh`). Never fight this by force-adding files
from those directories — if something there needs to be shared, it belongs in
`00_setup.sh`'s download/build logic, not committed as a binary blob.

## Pipeline execution order

`scripts/00` through `07` run in order (`01_verify` → `02_download` →
`03_dada2` → `04_kraken2_16s` → `05_build_tables` → `06_analysis` →
`07_compare_auroc`). `09_find_datasets.sbatch` and `10_inspect_candidates.sh` are
a separate, independent track (dataset discovery for phase 2) — they don't feed
into 01-07 automatically; a shortlisted accession has to be manually wired into
02/03/04/05/06 before it participates in the actual analysis. See "Adding a new
cohort" below.

`06_analysis.sbatch` is the one non-obvious piece: it's a generic runner over a
`MODULES` bash array (name | R script path | Singularity image | param | phase),
scaffolding each into a numbered `$PIPE/NN_Name/` directory and running the
original repo's R scripts unmodified inside the authors' own containers — **except**
for one patch applied at scaffold time, never to the reference clone: every bare
`fit<-anova(lme(...))` call gets wrapped in `tryCatch` (the authors already did
this themselves in `BS_kraken2_16S.R` for a known singular-convergence error; the
public Kraken2 substitute index surfaces the same failure in the other studies
too, so the patch generalizes their own fix rather than inventing a new one).
Phase-2 modules (`*Metagenomics`, `*pathway*`, `OppPathogens`) are listed but
skipped unless `RUN_PHASE2=1` — this repo has no metagenomic phase-2 work yet.

## Adding a new cohort (phase 2) — read this before touching 02-06

This is more than "write one new R script." Confirmed by reading the actual code,
not assumed:

1. **`02_download.sbatch` and `03_dada2.sbatch`** are SLURM-array `case`
   statements hardcoded to exactly the 4 phase-1 studies (`array=1-4`). A new
   cohort needs its own `case` branch in both (run list source, expected count,
   `truncLen`/primer/`trimLeft` params for DADA2).
2. **Per-cohort R scripts** are one hand-written file per study, not templated.
   The clone has **two parallel copies — do not confuse them**:
   - `Analysis/RScripts/{16S_TaxaClassification,DADA2_16S,kraken2_16S}/` — the
     **real, BioLockJ-style scripts**, the ones `06_analysis.sbatch`'s
     `MODULES` array actually references and runs (`pipeRoot <-
     dirname(dirname(getwd()))`-style path resolution). New cohort scripts
     must live here (well — in this repo's own tree once phase 2's scripts
     exist, mirroring this layout, then referenced by path from `MODULES`).
   - `Rcode/RYGB_IntegratedAnalysis/` — an **older/unused legacy copy** with
     the same filenames and near-identical logic, never invoked by anything in
     `scripts/`. Templating off this copy by mistake produces a script that
     looks right but is never executed. If you're reading a per-cohort script
     for reference, confirm the path starts with `Analysis/RScripts/` first.
   Pick the closest existing template (from `Analysis/RScripts/`) by
   timepoint shape, don't write from scratch:
   - **Ilhan's template** (3-level `Baseline`/`6M`/`12M` factor, `env_material`
     filter) — use for cohorts with a clean baseline + 2 follow-ups.
   - **Afshar's template** (binary `Pre`/`Post` factor, parses `ID`/`time` from
     a combined title string) — use for cohorts with one pre/post pair per
     subject.
3. **Four shared downstream scripts hardcode all 4 study names directly in
   code** (confirmed via `grep -rl "Afshar" Analysis/RScripts`):
   `combineCountTables.R`, `compareStudies_16S.R`, `compareStudies_SV.R`,
   `Heatmap.R`. Each needs a new cohort's variables/paths added by hand
   (literal `myT.<Study>` variables, matrix row-index math, `meta_all`
   construction in `combineCountTables.R` specifically). This is the part that
   actually breaks the "run the original repo's modules unmodified" framing —
   expected once you're past the original 4 cohorts, but be explicit about it
   in commit messages and README updates, don't silently absorb it into "ran
   the pipeline."
4. **New per-cohort metadata** (`cohorts/<Study>/metaData.txt` in this repo,
   not `$DATA_ROOT/RYGB_IntegratedAnalysis2020/input/...` — that tree is the
   *original authors'* data; ours doesn't belong there) needs a `Group`/`ID`
   (or `time`/`ID`) column pair derived from SRA `sample_alias`/`sample_title`,
   matching whichever template's expected column names. Verify readiness first
   with `scripts/10_inspect_candidates.sh <ACCESSION>` — checks read length →
   `truncLen`, specimen-type mixing (biopsy vs. fecal — filter to fecal only,
   as done for the ileocecal-resection cohort), and whether subject/timepoint
   parses from SRA alone or needs the source paper's own methods section.

**Timepoints never need to match across cohorts, including the new ones.** The
4 phase-1 cohorts already use 4 different scales (BS: 0/1/6 months; Afshar:
Pre/Post; Assal: Pre/3M/1Y/2Y; Ilhan: Baseline/6M/12M) — each is compared to its
own baseline via a per-study `Group`/`time` factor, never a common time axis. The
only real requirement is a recoverable pre-surgical baseline plus ≥1
post-surgical timepoint, repeat-sampled per subject.

**One study, one procedure.** If an accession bundles two different surgeries
under one BioProject (confirmed case: PRJNA1480144, an SDT arm and an ileostomy
arm under one accession), split it into two cohorts at integration time, each
following the closest single-procedure template — don't merge two procedures
into one `Study` label, since that blurs exactly the surgery-type signal this
whole project is testing for.

## SV/ASV-level comparison is region-locked — genus-level isn't

`compareStudies_SV.R` (exact-sequence/ASV-level cross-study comparison) only ever
compares **BS and Assal**, never Ilhan or Afshar, even though all 4 phase-1
studies have SV-level `MixedLinearModelResults` files. Why: BS and Assal both
use the 515F primer with `truncLen=200`, producing directly comparable 200bp
ASVs over the same amplified region; Ilhan (`truncLen=150`, no primer trim) and
Afshar (`trimLeft=20`, `truncLen=200`) produce different-length/region ASVs
that can't be meaningfully intersected at the exact-sequence level.
`compareStudies_16S.R` (genus-level) works across all 4 because genus
classification is region-agnostic — SV/ASV identity is not.

**This means a phase-2 cohort cannot default into SV comparison.** Verified
(2026-09-22) against real downloaded reads for all 4 phase-2 cohorts via the
515F-anchor check `03_dada2.sbatch` itself runs: none show a clean,
fixed-position anchor. IleocecalResection has zero matches across 5000
reads; Cholecystectomy/Ileostomy/SDT show scattered matches at varying
positions (159-190bp), the signature of coincidental short-motif hits within
the 16S sequence itself, not a real primer to strip. This confirms
`PRIMER=NONE`/`TRIMLEFT=0` is the *correct* setting for all 4 (no change
needed), and independently confirms none of them share BS/Assal's exact V4
protocol. The rewritten `analysis/RScripts/compareStudies/compareStudies_SV.R`
in this repo defaults to BS+Assal only (matching the original exactly) for
this reason; only add a phase-2 cohort to its `STUDIES` list after
confirming its real primer/region matches, not just because its `truncLen`
happens to also be 200.

## 16S vs. shotgun — don't assume from `scientific_name`

`AMPLICON` library strategy alone doesn't guarantee 16S (could be ITS, etc.), and
ENA's `scientific_name`/description sometimes reads "metagenome" for a 16S
amplicon sample (describes the sample source, not the sequencing method). Confirm
via the BioProject's actual description (`curl .../ena/browser/api/xml/<ACC>`,
look for "16S rDNA"/"16S rRNA" explicitly) or, failing that, the technical
fingerprint: `library_source=METAGENOMIC` + `library_selection=PCR` +
`AMPLICON` + a MiSeq 2×300bp or comparable read length is the standard 16S
V3-V4 signature.

## Python

Every Python entry point in this repo is a stdlib-only heredoc embedded directly
in `scripts/*.sbatch`/`scripts/*.sh` (`09_find_datasets.sbatch`,
`10_inspect_candidates.sh`, the convergence-guard patcher in `06_analysis.sbatch`)
— no third-party packages, by design, so a script stays self-contained and
portable across the cluster without a companion file to keep in sync. Don't
"fix" this by extracting them into a `src/` package unless there's an actual
reason (e.g., real unit tests are wanted) — `requirements.txt`/`pyproject.toml`
intentionally carry no runtime dependencies to match. If dataset-discovery logic
ever grows real test coverage, that's when `[tool.pytest.ini_options]` /
`[tool.coverage.run]` sections earn a place in `pyproject.toml` — not before.

## Verification philosophy

There's no pytest suite, and none is planned for the pipeline scripts themselves
— the actual verification loop is `07_compare_auroc.sbatch` and the comparison
block at the end of `06_analysis.sbatch`: rebuild a result from raw reads,
diff it directly against the original repo's own published output, and report
the delta rather than asserting a hardcoded expected value. Any new cohort
integration should follow the same pattern where a ground truth exists (the
source paper's own reported numbers), not invent unit tests for R modules that
belong to the original authors.

## ECC toolkit

A subset of the `ecc` agent/skill toolkit (cloned separately at
`/users/rgupta25/ecc`, sibling to this repo) is manually copied into
`~/.claude/agents/` (`python-reviewer`) and `~/.claude/skills/`
(`tdd-workflow`, `search-first`, `security-review`) — not a full plugin
install, since most of ECC's 68 agents (TypeScript/Kotlin/Java-focused) don't
apply to an R/Python/SLURM bioinformatics pipeline. If more of it is wanted
later, copy specific `agents/*.md` or `skills/<name>/` directories the same way
rather than running the full `ecc@ecc` plugin install.
