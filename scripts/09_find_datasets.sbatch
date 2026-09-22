#!/bin/bash
#SBATCH --job-name=rygb_scan
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/%u/rygb/logs/scan_%j.out
#SBATCH --error=/scratch/%u/rygb/logs/scan_%j.out
#
# Find GI-surgery 16S datasets that can be DROPPED INTO this pipeline.
#
# A dataset is only usable if all of these hold, so each is checked automatically:
#   1. AMPLICON + Illumina        — PacBio/454/IonTorrent need a separate DADA2 branch
#   2. read length >= 150 bp      — truncLen must reach 150 (Ilhan) or 200 (the others)
#   3. human gut                  — not mouse, not mucosa-only, not environmental
#   4. repeat sampling per subject — the mixed models need patient as a random effect
#   5. subject AND timepoint recoverable from SRA metadata alone
#                                 — no supplementary-table archaeology
#   6. a pre-surgical baseline plus >= 1 post-surgical timepoint
#                                 — each timepoint is compared to its own baseline,
#                                   so matching 1M/3M exactly is not required
#   7. GI-surgery context in the study title
#
# Output: ranked shortlist at $OUT/shortlist.tsv plus per-study run tables.

set -o pipefail
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
OUT=${OUT:-$HOME/rygb/dataset_scan}
mkdir -p "$OUT" "$DATA_ROOT/logs"
module load python/3.12.3

python3 - "$OUT" <<'PYEOF'
import json, re, sys, time, urllib.parse, urllib.request, os, collections

OUT = sys.argv[1]
API = "https://www.ebi.ac.uk/ena/portal/api"
FIELDS = ("run_accession,sample_accession,sample_alias,sample_title,library_name,"
          "library_strategy,library_layout,instrument_model,read_count,base_count,scientific_name")

# surgery vocabulary — gut / GI procedures, not bariatric only
KEYWORDS = [
    "gastrectomy", "gastric bypass", "sleeve gastrectomy", "bariatric surgery",
    "bariatric", "colectomy", "colorectal surgery", "colorectal resection",
    "rectal resection", "colorectal cancer surgery", "ileostomy", "stoma",
    "ileal pouch", "pouchitis", "ileocecal resection", "ileocolonic resection",
    "bowel resection", "intestinal resection", "esophagectomy",
    "pancreaticoduodenectomy", "cholecystectomy", "appendectomy",
    "gastrointestinal surgery", "abdominal surgery", "after surgery",
    "postoperative", "post-surgery", "surgical resection",
]

# accessions already known from the literature — checked first
SEEDS = {
    "PRJEB90497": "MedComm 2026 - colorectal surgery, pre-op + POD5-6, V3-V4",
}

NON_HUMAN = re.compile(r"\b(mouse|mice|murine|rat|porcine|pig|swine|chicken|bovine|"
                       r"canine|zebrafish|soil|marine|sediment|wastewater|rumen)\b", re.I)
B = r"(?<![A-Za-z0-9])"; E = r"(?![A-Za-z0-9])"
BASELINE = re.compile(B + r"(pre|pre-?op(erative)?|before|baseline|base|t0|d0|day ?0|v0|visit ?0)" + E, re.I)
TIMEPT = re.compile(B + r"(pre-?op(erative)?|post-?op(erative)?|baseline|before|after|"
                    r"pod ?\d+|[tvm]\d{1,2}|\d{1,3} ?(d|day|days|w|wk|week|weeks|"
                    r"m|mo|month|months|y|yr|year|years))" + E, re.I)

def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            if i == tries - 1:
                return f"__ERROR__ {e}"
            time.sleep(3 * (i + 1))

def search_studies(kw):
    q = urllib.parse.quote(f'study_title="*{kw}*"')
    txt = get(f"{API}/search?result=study&query={q}"
              f"&fields=study_accession,study_title&format=tsv&limit=300")
    if txt.startswith("__ERROR__") or not txt.strip():
        return []
    rows = []
    for line in txt.strip().split("\n")[1:]:
        p = line.split("\t")
        if len(p) >= 2:
            rows.append((p[0], p[1]))
    return rows

def runs_for(acc):
    path = os.path.join(OUT, f"{acc}.tsv")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        txt = open(path).read()
    else:
        txt = get(f"{API}/filereport?accession={acc}&result=read_run&fields={FIELDS}&format=tsv")
        if txt.startswith("__ERROR__"):
            return []
        open(path, "w").write(txt)
    lines = txt.strip().split("\n")
    if len(lines) < 2:
        return []
    hdr = lines[0].split("\t")
    return [dict(zip(hdr, l.split("\t"))) for l in lines[1:] if l.strip()]

def label_of(r):
    for k in ("sample_title", "sample_alias", "library_name"):
        v = (r.get(k) or "").strip()
        if v and v.lower() not in ("na", "none", "-"):
            return v
    return ""

def subject_timepoint(label):
    """split a sample label into (subject, timepoint) — the Palleja 'E02_1Y' pattern"""
    m = TIMEPT.search(label)
    if not m:
        return None, None
    tp = m.group(0).lower()
    subj = (label[:m.start()] + label[m.end():])
    subj = re.sub(r"[^A-Za-z0-9]+", "", subj).lower()
    return (subj or None), tp

def assess(acc, title, rows):
    n = len(rows)
    strat = {r.get("library_strategy", "") for r in rows}
    inst = {r.get("instrument_model", "") for r in rows}
    sci = {(r.get("scientific_name") or "").lower() for r in rows}
    reads = [int(r["read_count"]) for r in rows if (r.get("read_count") or "").isdigit()]
    bases = [int(r["base_count"]) for r in rows if (r.get("base_count") or "").isdigit()]
    paired = any((r.get("library_layout") or "") == "PAIRED" for r in rows)
    spot = (sum(bases) / sum(reads)) if reads and sum(reads) else 0
    readlen = spot / 2 if paired else spot

    labels = [label_of(r) for r in rows]
    st = [subject_timepoint(l) for l in labels]
    subjects = collections.Counter(s for s, _ in st if s)
    tps = collections.Counter(t for _, t in st if t)
    per_subj = (sum(subjects.values()) / len(subjects)) if subjects else 0
    has_base = any(BASELINE.search(l) for l in labels)

    reasons, score = [], 0
    if "AMPLICON" in strat and any("illumina" in i.lower() or "miseq" in i.lower()
                                   or "nextseq" in i.lower() or "novaseq" in i.lower()
                                   or "hiseq" in i.lower() for i in inst):
        score += 25
    else:
        reasons.append(f"not Illumina amplicon ({'/'.join(sorted(strat))}; {'/'.join(sorted(inst))})")
    if readlen >= 200: score += 10
    elif readlen >= 150: score += 7
    elif readlen: reasons.append(f"reads ~{readlen:.0f}bp, too short for truncLen 150")
    if n >= 40: score += 15
    elif n >= 20: score += 10
    else: reasons.append(f"only {n} runs")
    if len(subjects) >= 10 and per_subj >= 2: score += 20
    elif per_subj >= 2: score += 10
    else: reasons.append("no repeat sampling detectable from metadata")
    if has_base and len(tps) >= 2: score += 20
    elif len(tps) >= 2: score += 8; reasons.append("timepoints present but no clear baseline token")
    else: reasons.append("subject/timepoint not encoded in sample titles")
    if any("gut" in s or "feca" in s or "faeca" in s or "stool" in s for s in sci): score += 10
    elif any("human" in s for s in sci): score += 5
    if NON_HUMAN.search(title): score = 0; reasons.append("non-human")

    return dict(acc=acc, title=title, score=score, runs=n,
                strategy="/".join(sorted(strat)), instrument="/".join(sorted(inst))[:40],
                readlen=round(readlen), subjects=len(subjects),
                per_subject=round(per_subj, 1),
                timepoints=",".join(sorted(tps)[:8]), baseline=has_base,
                reasons="; ".join(reasons))

# ------------------------------------------------------------------ run
seen, results = set(), []

print("=" * 70)
print("A. Seed accessions from the literature")
print("=" * 70)
for acc, note in SEEDS.items():
    rows = runs_for(acc)
    if not rows:
        print(f"  {acc}: no runs returned"); continue
    a = assess(acc, note, rows); results.append(a); seen.add(acc)
    print(f"  {acc} score {a['score']:3d} | {a['runs']} runs | {a['strategy']} | "
          f"~{a['readlen']}bp | {a['subjects']} subjects x {a['per_subject']}")
    if a["reasons"]: print(f"      blockers: {a['reasons']}")

print()
print("=" * 70)
print("B. Keyword scan of ENA study titles")
print("=" * 70)
for kw in KEYWORDS:
    hits = search_studies(kw)
    print(f"\n--- '{kw}': {len(hits)} studies")
    for acc, title in hits:
        if acc in seen: continue
        seen.add(acc)
        if NON_HUMAN.search(title): continue
        rows = runs_for(acc)
        if len(rows) < 12: continue
        a = assess(acc, title, rows)
        results.append(a)
        if a["score"] >= 60:
            print(f"   ** {acc} score {a['score']:3d} | {a['runs']} runs | ~{a['readlen']}bp | "
                  f"{a['subjects']} subj x {a['per_subject']} | {title[:60]}")
            print(f"      timepoints: {a['timepoints']}")
        time.sleep(0.25)

# ------------------------------------------------------------------ report
results.sort(key=lambda d: -d["score"])
cols = ["score","acc","runs","subjects","per_subject","readlen","strategy",
        "instrument","baseline","timepoints","title","reasons"]
with open(os.path.join(OUT, "shortlist.tsv"), "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in results:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")

print()
print("=" * 70)
print("TOP CANDIDATES  (score >= 70 = drop-in; 50-69 = usable with work)")
print("=" * 70)
for r in results[:25]:
    if r["score"] < 40: break
    print(f"\n{r['score']:3d}  {r['acc']}  {r['title'][:70]}")
    print(f"     {r['runs']} runs | {r['subjects']} subjects x {r['per_subject']} samples "
          f"| ~{r['readlen']}bp | {r['instrument']}")
    print(f"     timepoints seen: {r['timepoints']}  (baseline token: {r['baseline']})")
    if r["reasons"]: print(f"     caveats: {r['reasons']}")

print(f"\nfull table: {os.path.join(OUT,'shortlist.tsv')}")
print("""
Still needs a human for each survivor:
  - 16S region and primers  -> decides cutadapt vs trimLeft, from the paper's methods
  - read length             -> decides truncLen (200 like BS/Assal/Afshar, 150 like Ilhan)
  - perioperative antibiotics -> a confounder to record, not a disqualifier
  - the surgery itself      -> a Roux-en-Y reconstruction after gastrectomy tests the
                               same rewiring as RYGB in different patients; a colectomy
                               is the negative control for "any gut surgery does this"
""")
PYEOF
