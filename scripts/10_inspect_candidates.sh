#!/bin/bash
# 10_inspect_candidates.sh — decide which scanned studies are actually usable.
# Reads the run tables cached by 09_find_datasets.sbatch, so it re-downloads nothing.
# Runs in seconds on a login node:
#     bash 10_inspect_candidates.sh 2>&1 | tee ~/rygb/candidate_inspect.txt
#
# Does two things:
#   1. RE-SCORES every cached study with the non-human filter applied to scientific_name
#      (the scan's title-only filter let shrimp and pika stomachs through), and stops
#      treating "subject/timepoint missing from SRA" as disqualifying — PRJNA668472,
#      the BS cohort already in this pipeline, fails that check, because its metadata
#      was built from the paper. It is now reported as "metadata needed from paper".
#   2. DEEP-INSPECTS named accessions: every distinct sample label, the subject x
#      timepoint layout, read lengths and layout — everything needed to decide whether
#      a dataset drops in, and with what truncLen.

set -o pipefail
OUT=${OUT:-$HOME/rygb/dataset_scan}
module load python/3.12.3 2>/dev/null

python3 - "$OUT" "$@" <<'PYEOF'
import os, re, sys, glob, collections

OUT = sys.argv[1]
FOCUS = sys.argv[2:] or [
    "PRJNA1188648",   # cholecystectomy, baseline/6M/12M
    "PRJEB39382",     # bariatric, day-level timepoints, 258 runs
    "PRJNA514452",    # ileocecal resection, 509 runs
    "PRJNA1130449",   # ileostomy closure + probiotic
    "PRJNA655569",    # appendectomy
    "PRJNA727576",    # bariatric, 358 runs
    "PRJNA1480144",   # SDT + ileostomy
    "PRJNA703442",    # CRC prolonged ileus
    "PRJNA668472",    # the BS cohort already in the pipeline — the control
]

B = r"(?<![A-Za-z0-9])"; E = r"(?![A-Za-z0-9])"
TIMEPT = re.compile(B + r"(pre-?op(erative)?|post-?op(erative)?|baseline|before|after|"
                    r"pod ?\d+|[tvm]\d{1,2}|\d{1,3} ?(d|day|days|w|wk|week|weeks|"
                    r"m|mo|month|months|y|yr|year|years))" + E, re.I)
BASELINE = re.compile(B + r"(pre|pre-?op(erative)?|before|baseline|base|t0|d0|v0)" + E, re.I)
NONHUMAN_SCI = re.compile(r"(shrimp|pika|mouse|mus |rat |murine|pig|swine|bovine|chicken|"
                          r"fish|aquacultur|soil|marine|sediment|water|plant|insect|"
                          r"penaeus|marsupenaeus)", re.I)

def load(acc):
    p = os.path.join(OUT, f"{acc}.tsv")
    if not os.path.exists(p): return None
    lines = open(p).read().strip().split("\n")
    if len(lines) < 2: return None
    hdr = lines[0].split("\t")
    return [dict(zip(hdr, l.split("\t"))) for l in lines[1:] if l.strip()]

LABEL_FIELDS = ("sample_title", "sample_alias", "library_name", "sample_accession")

def field_cardinality(rows):
    """how many distinct values each candidate label field has"""
    out = {}
    for k in LABEL_FIELDS:
        vals = [(r.get(k) or "").strip() for r in rows]
        vals = [v for v in vals if v and v.lower() not in ("na", "none", "-")]
        out[k] = (len(set(vals)), vals[0] if vals else "")
    return out

def best_field(rows):
    """pick the most informative field, not merely the first non-empty one.
       PRJNA668472 (the BS cohort) carries boilerplate in sample_title and the real
       IDs in sample_alias — taking sample_title first made it look unusable."""
    card = field_cardinality(rows)
    # prefer a field that both varies and parses into subject+timepoint
    scored = []
    for k, (n, _) in card.items():
        if n <= 1: continue
        vals = [(r.get(k) or "").strip() for r in rows]
        parsed = sum(1 for v in vals if TIMEPT.search(v))
        scored.append((parsed, n, k))
    if not scored:
        return None
    scored.sort(reverse=True)
    return scored[0][2]

def label(r, field=None):
    if field:
        return (r.get(field) or "").strip()
    for k in LABEL_FIELDS:
        v = (r.get(k) or "").strip()
        if v and v.lower() not in ("na", "none", "-"): return v
    return ""

def split_label(l):
    m = TIMEPT.search(l)
    if not m: return None, None
    return (re.sub(r"[^A-Za-z0-9]+", "", l[:m.start()] + l[m.end():]).lower() or None,
            m.group(0).lower())

def readlen(rows):
    reads = sum(int(r["read_count"]) for r in rows if (r.get("read_count") or "").isdigit())
    bases = sum(int(r["base_count"]) for r in rows if (r.get("base_count") or "").isdigit())
    paired = any(r.get("library_layout") == "PAIRED" for r in rows)
    if not reads: return 0, paired
    return (bases / reads) / (2 if paired else 1), paired

print("=" * 72)
print("1. RE-SCORED — non-human removed by scientific_name, metadata-gap downgraded")
print("=" * 72)
rows_out = []
for p in sorted(glob.glob(os.path.join(OUT, "PRJ*.tsv"))):
    acc = os.path.basename(p)[:-4]
    rows = load(acc)
    if not rows: continue
    sci = collections.Counter((r.get("scientific_name") or "?") for r in rows)
    top_sci = sci.most_common(1)[0][0]
    if NONHUMAN_SCI.search(" ".join(sci)):
        continue                                     # shrimp, pika, mouse...
    strat = {r.get("library_strategy") for r in rows}
    if "AMPLICON" not in strat: continue
    rl, paired = readlen(rows)
    if rl < 140: continue                            # cannot reach truncLen 150
    bf = best_field(rows)
    labels = [label(r, bf) for r in rows]
    st = [split_label(l) for l in labels]
    subj = collections.Counter(s for s, _ in st if s)
    tps = sorted({t for _, t in st if t})
    per = (sum(subj.values()) / len(subj)) if subj else 0
    meta_ok = len(subj) >= 5 and per >= 1.8
    rows_out.append((len(rows), acc, top_sci, round(rl), len(subj), round(per, 1),
                     meta_ok, ",".join(tps[:6])))
rows_out.sort(key=lambda r: (-r[6], -r[0]))
print(f"{'runs':>5} {'accession':14} {'source':26} {'bp':>4} {'subj':>5} {'/subj':>5} "
      f"{'metadata':10} timepoints")
for n, acc, sci, rl, ns, per, ok, tps in rows_out:
    print(f"{n:5d} {acc:14} {sci[:26]:26} {rl:4d} {ns:5d} {per:5.1f} "
          f"{'in SRA' if ok else 'from paper':10} {tps}")

print()
print("=" * 72)
print("2. DEEP INSPECTION")
print("=" * 72)
for acc in FOCUS:
    rows = load(acc)
    print(f"\n{'-'*72}\n{acc}")
    if not rows:
        print("   not in cache — add it to 09_find_datasets.sbatch SEEDS and rerun")
        continue
    sci = collections.Counter((r.get("scientific_name") or "?") for r in rows)
    inst = collections.Counter((r.get("instrument_model") or "?") for r in rows)
    lay = collections.Counter((r.get("library_layout") or "?") for r in rows)
    rl, paired = readlen(rows)
    print(f"   {len(rows)} runs | {dict(sci)} | {dict(inst)} | {dict(lay)}")
    print(f"   mean read length ~{rl:.0f} bp  ->  truncLen "
          f"{'200 (like BS/Assal/Afshar)' if rl >= 210 else '150 (like Ilhan)' if rl >= 155 else 'TOO SHORT'}")
    card = field_cardinality(rows)
    print("   label fields (distinct values | example):")
    for k, (n, ex) in card.items():
        print(f"      {k:18} {n:5d} | {ex[:60]}")
    bf = best_field(rows)
    print(f"   -> using: {bf or 'NONE — no field parses into subject+timepoint'}")
    labels = [label(r, bf) for r in rows]
    uniq = sorted(set(labels))
    print(f"   {len(uniq)} distinct labels; first 24:")
    for l in uniq[:24]:
        print(f"      {l}")
    st = [split_label(l) for l in labels]
    subj = collections.Counter(s for s, _ in st if s)
    tps = collections.Counter(t for _, t in st if t)
    if subj:
        print(f"   parsed: {len(subj)} subjects, {len(tps)} timepoints")
        print(f"   timepoint counts: {dict(tps.most_common(12))}")
        multi = sum(1 for v in subj.values() if v >= 2)
        print(f"   subjects with >=2 samples: {multi}/{len(subj)}")
        print(f"   baseline token present: {any(BASELINE.search(l) for l in labels)}")
    else:
        print("   no subject/timepoint parsed -> metadata must come from the paper's")
        print("   supplementary table, exactly as the authors did for the BS cohort")
PYEOF
