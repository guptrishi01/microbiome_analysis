# COCKTAIL_GBP cohort — metadata gap (more severe than Penney/Oseberg)

`metaData.txt` has REAL `Run`, `BioProject`, and `sample_alias` for all 408
runs under PRJNA904651 (confirmed against the live ENA run report), but
`time`, `ID`, and `Arm` are ALL `UNKNOWN` — not placeholders I chose not to
fill in, but genuinely not exposed anywhere in public sequence-archive
metadata. `sample_alias` is an anonymized hash (e.g. `b1951117`), not a
patient-ID/timepoint string like Penney's or Oseberg's, and `sample_title`
is blank for every run.

This isn't an oversight by the authors — it's deliberate. Gopalakrishnan et
al. 2025's own Data Availability Statement (Mucosal Immunology 18(4):825-835,
doi:10.1016/j.mucimm.2025.04.002) says: "De-identified individual data and
metadata... may be shared upon reasonable request evaluated in the COCKTAIL
Steering Committee." The subject/timepoint/arm mapping was intentionally
kept out of the public archive.

**Before this cohort (GBP arm only — the VLED diet arm is explicitly excluded
per this project's scope decision) can be processed at all:**
1. Either request the de-identified metadata from the COCKTAIL Steering
   Committee via the corresponding author (Jøran Hjelmesæth,
   joran.hjelmeseth@siv.no, per the paper), or
2. Check whether the paper's own Supplementary Tables (not the sequence
   archive) include a run-accession-to-sample mapping — unconfirmed either
   way; I have not been able to check this directly, so I'm flagging it as a
   possibility to check, not a tested alternative.
3. Once a Run -> {time, ID, Arm} mapping exists from either source, fill in
   this file and remove the UNKNOWN placeholders.

Until then, this cohort cannot be meaningfully taken past raw-read download —
`03b_kneaddata.sbatch` could run for QC purposes, but
`04b_kraken2_wgs.sbatch`'s output would have no way to be grouped by
timepoint, and the per-cohort MLM script cannot be written without knowing
which samples are pre/post and which patient each belongs to.

**Re-checked 2026-11-14** (no new per-run resolution found): a web search
surfaced the paper's ScienceDirect listing
(https://www.sciencedirect.com/science/article/pii/S1933021925000406), whose
abstract gives one new, usable-as-context (not actionable) public fact: the
trial is "non-randomized, single-center... 88 participants with severe
obesity, split between gastric bypass and a very low energy diet [VLED] that
was matched for weight loss over 6 weeks," followed 2 years. That's a total
n=88 across BOTH arms (not per-arm, not per-run) — it does not resolve any
Run's time/ID/Arm and does not change anything in metaData.txt.

I could not re-fetch the article's actual Data Availability Statement text
directly in this pass (no page-fetch tool was available, only web search,
which returns snippets/summaries, not full article text) — the DAS quote and
contact above are carried over from earlier session research, not
independently re-verified this pass. Before relying on the
joran.hjelmeseth@siv.no contact path, have a pass with a working fetch tool
(or WebFetch) re-read https://doi.org/10.1016/j.mucimm.2025.04.002 directly
to confirm that quote verbatim.

**INDEPENDENTLY RE-VERIFIED 2026-11-14**, directly from the primary source:
the user supplied the full PDF (`1-s2.0-S1933021925000406-main.pdf`, read and
then deleted per this repo's PDF-handling convention — never committed).
Its "Availability of Data and Materials" section reads, verbatim: "Fecal
whole metagenomic sequencing data will be made available through the NCBI
sequence read archive PRJNA904651. De-identified individual data and
metadata, and processed metabolomic data may be shared upon reasonable
request evaluated in the COCKTAIL Steering Committee. Requested analyses
must be within the boundaries of the patient informed consent and approved
by the Regional Committee for Medical and Health Research Ethics (REK
sørøst A, Norway). Applications to the ethical committee will be submitted
in cooperation with the principal investigator Professor Jøran Hjelmesæth
(joran.hjelmeseth@siv.no)." This matches the earlier-session quote exactly
— confirmed, not carried over unverified anymore. The paper's own Tables
S1-S4 (pathway/metabolite associations) were also confirmed to contain no
per-patient ID/arm/timepoint mapping. COCKTAIL_GBP remains genuinely
blocked on the Steering Committee request; there is no other public route.
