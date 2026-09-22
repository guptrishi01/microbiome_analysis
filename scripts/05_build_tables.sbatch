#!/bin/bash
#SBATCH --job-name=rygb_tables
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/%u/rygb/logs/tables_%j.out
#SBATCH --error=/scratch/%u/rygb/logs/tables_%j.out
#
# Stage 5: mpa reports -> *_taxaCount_norm_Log10_<Rank>.tsv, in the repo's layout and filenames.
#   - keeps only d__Bacteria and d__Archaea (the authors' custom DB held nothing else),
#     dropping Eukaryota/Viruses BEFORE per-sample totals are computed
#   - normalization: log10( count / sample_total * mean(sample_total) + 1 )
#   - copies the hand-curated metadata from the reference clone
#   - compares every regenerated table against the repo's own

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
REPO=$DATA_ROOT/RYGB_IntegratedAnalysis2020
COHORTS=$SW_ROOT/cohorts
module load python/3.12.3

cat > "$DATA_ROOT/work/build_tables.py" <<'PYEOF'
import os, sys, glob, math, shutil, statistics

DATA_ROOT = sys.argv[1]
REPO      = sys.argv[2]
COHORTS   = sys.argv[3]

RANKS = [('p','Phylum'), ('c','Class'), ('o','Order'), ('f','Family'), ('g','Genus'), ('s','Species')]
KEEP_DOMAINS = {'d__Bacteria', 'd__Archaea'}

STUDIES = {
    'BS16S':  dict(repodir='RYGB_BS',         prefix='RYGB_BS_2020Oct04',      n=61),
    'Assal':  dict(repodir='RYGB_Assal2020',  prefix='RYGB_Assal_2020Apr23',   n=71),
    'Afshar': dict(repodir='RYGB_Afshar2018', prefix='RYGB_Afshar_2020Apr23',  n=50),
    'Ilhan':  dict(repodir='RYGB_Ilhan2020',  prefix='RYGB_Ilhan_2020Apr23',   n=82),
    # Phase 2 — no entry in the original repo; metadata comes from cohorts/<key>/ in this
    # repo instead (see the metadata-copy fallback below), not $REPO.
    'Cholecystectomy':    dict(repodir='RYGB_Cholecystectomy',    prefix='RYGB_Cholecystectomy_2026',    n=30,  own_metadata=True),
    'IleocecalResection': dict(repodir='RYGB_IleocecalResection', prefix='RYGB_IleocecalResection_2026', n=189, own_metadata=True),
    'Ileostomy':          dict(repodir='RYGB_Ileostomy',          prefix='RYGB_Ileostomy_2026',          n=65,  own_metadata=True),
    'SDT':                dict(repodir='RYGB_SDT',                prefix='RYGB_SDT_2026',                n=45,  own_metadata=True),
}

def parse_mpa(path):
    """deepest-field rank -> {taxon: clade count}, bacteria and archaea only"""
    out = {r: {} for r, _ in RANKS}
    dropped = 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            try:
                lineage, count = line.rsplit('\t', 1)
                count = int(count)
            except ValueError:
                continue
            fields = lineage.split('|')
            if fields[0] not in KEEP_DOMAINS:
                dropped += count
                continue
            last = fields[-1]
            if len(last) < 4 or last[1:3] != '__':
                continue
            rank = last[0]
            if rank in out:
                name = last[3:].replace(chr(34),"").replace(chr(39),"").replace("#","")
                out[rank][name] = out[rank].get(name, 0) + count
    return out, dropped

def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return float('nan')
    mx, my = sum(xs)/n, sum(ys)/n
    sxy = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
    sxx = sum((x-mx)**2 for x in xs)
    syy = sum((y-my)**2 for y in ys)
    if sxx == 0 or syy == 0:
        return float('nan')
    return sxy / math.sqrt(sxx*syy)

def read_repo_table(path):
    """SAMPLE_ID table -> (taxa, {sample: {taxon: value}})"""
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        taxa = header[1:]
        rows = {}
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            rows[parts[0]] = {t: float(v) for t, v in zip(taxa, parts[1:])}
    return taxa, rows

for study, cfg in STUDIES.items():
    print(f"\n{'='*20} {study}")
    mpas = sorted(glob.glob(os.path.join(DATA_ROOT, 'work', 'kraken', study, '*.mpa')))
    print(f"mpa reports: {len(mpas)} (expect {cfg['n']})")
    if not mpas:
        print("  none found, skipping")
        continue

    samples = [os.path.basename(p)[:-4] for p in mpas]
    parsed, dropped = {}, 0
    for p, s in zip(mpas, samples):
        parsed[s], d = parse_mpa(p)
        dropped += d
    print(f"non-bacterial/archaeal reads dropped before normalization: {dropped:,}")

    outdir = os.path.join(DATA_ROOT, 'input', cfg['repodir'], 'Kraken2')
    os.makedirs(outdir, exist_ok=True)

    for rank, rank_name in RANKS:
        taxa = sorted({t for s in samples for t in parsed[s][rank]})
        totals = {s: sum(parsed[s][rank].values()) for s in samples}
        nonzero = [totals[s] for s in samples if totals[s] > 0]
        mean_total = statistics.mean(nonzero) if nonzero else 0.0

        outfile = os.path.join(outdir, f"{cfg['prefix']}_taxaCount_norm_Log10_{rank_name}.tsv")
        with open(outfile, 'w') as fh:
            fh.write('SAMPLE_ID\t' + '\t'.join(taxa) + '\n')
            for s in sorted(samples):
                tot = totals[s]
                vals = []
                for t in taxa:
                    c = parsed[s][rank].get(t, 0)
                    v = math.log10(c / tot * mean_total + 1) if tot > 0 and c > 0 else 0.0
                    vals.append(repr(v) if v else '0.0')
                fh.write(s + '\t' + '\t'.join(vals) + '\n')

        # ---- compare with the published table
        repofile = os.path.join(REPO, 'input', cfg['repodir'], 'Kraken2',
                                f"{cfg['prefix']}_taxaCount_norm_Log10_{rank_name}.tsv")
        msg = ''
        if os.path.exists(repofile):
            rtaxa, rrows = read_repo_table(repofile)
            shared_taxa = sorted(set(rtaxa) & set(taxa))
            _, orows = read_repo_table(outfile)
            shared_samples = sorted(set(rrows) & set(orows))
            xs, ys = [], []
            for s in shared_samples:
                for t in shared_taxa:
                    xs.append(rrows[s][t]); ys.append(orows[s][t])
            r = pearson(xs, ys)
            nz = [(x, y) for x, y in zip(xs, ys) if x > 0 or y > 0]
            rnz = pearson([x for x, _ in nz], [y for _, y in nz]) if nz else float('nan')
            msg = (f"| repo {len(rtaxa):5d} taxa, shared {len(shared_taxa):5d} "
                   f"| r={r:.3f} (nonzero-only r={rnz:.3f}, n={len(nz):,})")
        print(f"  {rank_name:8s} {len(taxa):5d} taxa  mean depth {mean_total:9.0f}  {msg}")

    # ---- metadata: the authors' own hand-curated files for phase-1 studies,
    # or our own cohorts/<study>/metaData.txt for phase-2 studies (own_metadata=True)
    copied = []
    if cfg.get('own_metadata'):
        src = os.path.join(COHORTS, study, 'metaData.txt')
        if os.path.exists(src):
            dst_dir = os.path.join(DATA_ROOT, 'input', cfg['repodir'], 'Metadata')
            os.makedirs(dst_dir, exist_ok=True)
            shutil.copy2(src, os.path.join(dst_dir, 'metaData.txt'))
            copied.append('Metadata/metaData.txt (from cohorts/)')
        else:
            print(f"  !! no metadata found at {src} — build it before running 06_analysis")
    else:
        for sub, pattern in (('Metadata', '*'), ('Kraken2', '*.txt')):
            src = os.path.join(REPO, 'input', cfg['repodir'], sub)
            if not os.path.isdir(src):
                continue
            dst = os.path.join(DATA_ROOT, 'input', cfg['repodir'], sub)
            os.makedirs(dst, exist_ok=True)
            for f in glob.glob(os.path.join(src, pattern)):
                if 'speciesMetadata' in os.path.basename(f):   # generated by the R stage
                    continue
                shutil.copy2(f, dst)
                copied.append(f"{sub}/{os.path.basename(f)}")
    print("  metadata copied:", ', '.join(copied) if copied else 'none')

print("\nDone.")
PYEOF

python3 "$DATA_ROOT/work/build_tables.py" "$DATA_ROOT" "$REPO" "$COHORTS"
rc=$?

echo
echo "=== regenerated input tree"
find "$DATA_ROOT/input" -type f | sort | sed 's/^/   /'
exit $rc
