# Penney_RYGB cohort — partial resolution, more patients may exist

`metaData.txt` here (25 rows, 10 patients: 012, 014, 015, 016, 021, 030, 039,
042, 046, 048) is a **confirmed-RYGB subset** of the full Penney et al. 2022
dataset (PRJNA473348, 120 runs / ~81 patients total). Resolved 2026-11-14 from
the paper's own Mendeley Data deposit (DOI 10.17632/t76nm3yfzh.3,
`Source Data/Figure 6/MixOmics_PairedRYGB.xlsx`'s `Sample_Name` column) — see
git history / session notes for the full verification (ID-space match against
ENA `sample_alias`, cross-check against the paper's own reported n=13 RYGB in
its "27 complete paired sets" figure).

**This is NOT the full RYGB arm.** The paper reports 13 complete-paired RYGB
patients; only 10 appear in this specific multi-omics source-data file (likely
because that 4-block DIABLO model required complete data across every omics
layer — serum, urine, faecal — a stricter filter than microbiome pairing
alone). The other ~3 RYGB patients, and the ~54 additional Penney patients
with no Arm information anywhere in the public Mendeley deposit (exhaustively
checked — see session notes), are simply not in this cohort. If you obtain a
fuller patient ID → Arm mapping later (from the corresponding authors,
n.penney@imperial.ac.uk or elaine.holmes@imperial.ac.uk), add those patients'
rows here (and to Penney_SG/ as appropriate) rather than starting over.
