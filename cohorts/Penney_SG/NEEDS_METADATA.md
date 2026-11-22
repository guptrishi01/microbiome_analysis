# Penney_SG cohort — partial resolution, more patients may exist

`metaData.txt` here (31 rows, 14 patients: 004, 005, 007, 011, 017, 018, 040,
041, 049, 051, 052, 056, 065, 093) is a **confirmed-SG (sleeve gastrectomy)
subset** of the full Penney et al. 2022 dataset (PRJNA473348, 120 runs / ~81
patients total). Resolved 2026-11-14 from the paper's own Mendeley Data
deposit (DOI 10.17632/t76nm3yfzh.3,
`Source Data/Figure 6/MixOmics_PairedVSG.xlsx`'s `Sample_Name` column) — see
git history / session notes for the full verification. The paper calls this
arm "VSG" (Vertical Sleeve Gastrectomy); mapped to this repo's standard "SG"
label for consistency with Oseberg_SG/Lazaro_SG/Chaudhari_SG.

**This matches the full SG arm exactly** — the paper reports 14
complete-paired VSG patients in its "27 complete paired sets" figure, and
this file has exactly 14. (Contrast with Penney_RYGB/, which is short 3
patients from the paper's reported 13 RYGB — see that cohort's own
NEEDS_METADATA.md.) The ~54 additional Penney patients with no Arm
information anywhere in the public Mendeley deposit (exhaustively checked —
see session notes) are not in this cohort either. If you obtain a fuller
patient ID → Arm mapping later (from the corresponding authors,
n.penney@imperial.ac.uk or elaine.holmes@imperial.ac.uk), add those patients'
rows here (and to Penney_RYGB/ as appropriate) rather than starting over.
