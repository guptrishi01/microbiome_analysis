# Hayashi_LSG cohort — metadata gap (time/ID resolved, Arm still open)

## Resolved (2026-11-14): `time` and `ID`

`metaData.txt`'s `time`/`ID` columns are no longer placeholders. The archive's
own exposed fields (`sample_alias`, `sample_title`) really do carry no
per-sample info (`sample_alias` is a generic `SampleNNN` string; `sample_title`
is the boilerplate MIxS description for every run) — but each BioSample's full
XML record (not surfaced by ENA's filereport, only via NCBI
`efetch?db=biosample`) carries a `source_material_id` attribute like `OBB8` or
`OBA19` that ENA drops. Decoded across all 38 BioSamples:

- The numeric suffix (1-16, 18-20 — 19 distinct values, patient 17 absent) is
  the patient ID. Every one of the 19 has exactly one `OBA*` and one `OBB*`
  sample — confirmed programmatically, no exceptions.
- The `A`/`B` letter is timepoint, not arm: for every one of the 19 patients,
  the `OBB*` sample's `collection_date` is chronologically BEFORE the
  matching `OBA*` sample's (gap 177-418 days, nominally ~6 months, with
  real-world clinical-followup variance) — 19/19 consistent, no exceptions.
  So `B` = `Pre` (baseline), `A` = `Post` (follow-up). This also lines up with
  the paper's own terminology: "OB" = Obese-group label (same prefix used for
  the "OB"/"HD" obese-vs-healthy-donor split in the paper's Supplementary
  Table S1), not a procedure code.
- Verified via: `curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=<SAMN...>&rettype=full&retmode=xml"`
  for all 38 BioSample accessions (joined to Run via the ENA filereport's
  `sample_accession` field) — each record's
  `<Attribute attribute_name="source_material_id">` field.

`metaData.txt` now has real `Pre`/`Post` in `time` and the patient number
(1-20, skipping 17) in `ID` for all 38 runs.

## Still UNKNOWN: `Arm` (LSG vs. IGBP)

The paper (Hayashi T, et al. "Escherichia coli as a gut microbial marker of
obesity and its reduction following bariatric treatment." J Gastroenterol.
2026;61(6):741-749. doi:10.1007/s00535-026-02362-z — open access) confirms 19
patients (13 LSG, 6 IGBP non-surgical — excluded from this cohort per this
project's scope decision), sampled at baseline and 6 months, 38 samples total
— matching exactly. But neither the full open-access article text, nor any
of its 3 published tables (Table 1: aggregate obese-vs-healthy-donor stats;
Table 2: aggregate pre/post BMI, "Bariatric procedure (IGBP:LSG) 6:13" as a
ratio only; Table 3: aggregate clinical-outcome stats), nor its Supplementary
file (`535_2026_2362_MOESM1_ESM.docx` — contains only Supplementary Fig. 1
and Tables S1/S2, both KEGG-pathway/clinical-correlation coefficient tables,
no per-patient data) names which of the 19 patient numbers underwent LSG vs.
IGBP. The BioSample XML's `source_material_id` doesn't encode it either (it's
patient+timepoint, not procedure, per above). This is a genuine, confirmed
gap, not an unchecked possibility — all 3 named public sources were checked
directly.

**Before this cohort can be processed:**
1. Contact the corresponding author (Noriho Iida, niida@m-kanazawa.jp, per
   the paper) for the patient-ID-to-{LSG, IGBP} mapping. Patient IDs to ask
   about: 1-16, 18-20 (this cohort's own numbering, derived from each
   BioSample's `source_material_id` — give them this numbering or the
   Run/SAMN accessions so they can match it to their own records).
2. Once resolved: fill in real `Arm` values (filter to LSG-only — 13 patients
   x 2 timepoints = 26 of the 38 runs — before this cohort enters the
   pipeline, per this project's scope decision to exclude IGBP).

Until then, `Hayashi_LSG_kraken2_Metagenomics.R`'s `stop()` guard (checks
`time`/`ID`/`Arm` for `UNKNOWN`) will still correctly block this cohort —
`time`/`ID` no longer trip it, but `Arm` still does, which is now the only
thing blocking this cohort.
