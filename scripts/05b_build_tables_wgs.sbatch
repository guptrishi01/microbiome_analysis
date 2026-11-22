#!/bin/bash
#SBATCH --job-name=rygb_tables_wgs
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/%u/rygb/logs/tables_wgs_%j.out
#SBATCH --error=/scratch/%u/rygb/logs/tables_wgs_%j.out
#
# Stage 5 (WGS track): Bracken per-rank output -> <prefix>_taxaCount_norm_Log10_<Rank>.tsv,
# same filename convention and log10 normalization formula
# (log10(count/sample_total*mean(sample_total) + 1)) as:
#   - 05_build_tables.sbatch (16S track, this repo)
#   - the original authors' own RYGB_BS_metagenomics_2020Oct01_taxaCount_norm_Log10_*.tsv
#     (input/RYGB_BS_Metagenomics/ in the read-only clone) — so BS-WGS and these 4 new
#     cohorts are drop-in interchangeable inputs to the same per-cohort MLM script shape.
#
# Uses Bracken's FLAT .bracken file per rank (name/taxonomy_id/taxonomy_lvl/
# kraken_assigned_reads/added_reads/new_est_reads/fraction_total_reads) for counts, joined
# against each sample's .k2report (already produced by 04b for Bracken's own input) to
# recover domain membership per taxid — Bracken's flat file itself has no parent-domain
# column, but the k2report's rank_code/indentation does encode the full lineage, so no
# second Kraken2 run is needed. Filters to Bacteria/Archaea only, dropped BEFORE per-sample
# totals are computed, exactly matching 05_build_tables.sbatch's 16S KEEP_DOMAINS behavior
# (human reads are additionally already removed upstream by kneaddata in 03b, unlike 16S's
# input, so the main thing this filter drops here is viral/unclassified, not human).
#
# Copies this repo's own cohorts/<Study>/metaData.txt alongside the tables, same as
# 05_build_tables.sbatch does for own_metadata=True cohorts.

set -o pipefail
SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
COHORTS=$SW_ROOT/cohorts
module load python/3.12.3

cat > "$DATA_ROOT/work/build_tables_wgs.py" <<'PYEOF'
import os, sys, glob, math, shutil, statistics

DATA_ROOT = sys.argv[1]
COHORTS   = sys.argv[2]

RANKS = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

# key -> (repodir under $DATA_ROOT/input/, output filename prefix, expected run count)
# Every one of these is "own_metadata=True"-equivalent: metadata comes from this repo's
# cohorts/<key>/metaData.txt, never from the read-only original-authors clone (which has
# no entry for any of them). BS-WGS is handled separately — see 02b_download_wgs.sbatch's
# header comment for why it's not reprocessed here at all.
STUDIES = {
    # Penney split 2026-11-14 once a PARTIAL Arm resolution was found (Mendeley deposit
    # linked from Penney et al. 2022's own Methods) — was one unsplit 'Penney' entry
    # (n=120), now two single-procedure cohorts per CLAUDE.md's "one study, one
    # procedure" rule. Only 25+31=56 of the original 120 runs are covered — the rest
    # have no recoverable Arm (see cohorts/Penney_RYGB/ and Penney_SG/NEEDS_METADATA.md).
    'Penney_RYGB':  dict(repodir='RYGB_Penney_RYGB_WGS',  prefix='RYGB_Penney_RYGB_2026_metagenomics',  n=25),
    'Penney_SG':    dict(repodir='RYGB_Penney_SG_WGS',    prefix='RYGB_Penney_SG_2026_metagenomics',    n=31),
    # Oseberg split 2026-11-14 once Arm was resolved (Zenodo source-data deposit linked
    # from Olsson et al. 2026's Data Availability Statement) — was one unsplit 'Oseberg'
    # entry (n=154), now two single-procedure cohorts per CLAUDE.md's "one study, one
    # procedure" rule.
    'Oseberg_RYGB': dict(repodir='RYGB_Oseberg_RYGB_WGS', prefix='RYGB_Oseberg_RYGB_2026_metagenomics', n=78),
    'Oseberg_SG':   dict(repodir='RYGB_Oseberg_SG_WGS',   prefix='RYGB_Oseberg_SG_2026_metagenomics',   n=76),
    'COCKTAIL_GBP': dict(repodir='RYGB_COCKTAIL_GBP_WGS',  prefix='RYGB_COCKTAIL_GBP_2026_metagenomics', n=408),
    'Hayashi_LSG':  dict(repodir='RYGB_Hayashi_LSG_WGS',   prefix='RYGB_Hayashi_LSG_2026_metagenomics',  n=38),
}

KEEP_DOMAINS = {'Bacteria', 'Archaea'}

def parse_k2report_domains(path):
    """k2report (standard, indented) -> {taxid: domain or None}, domain in KEEP_DOMAINS or None"""
    domain_by_taxid = {}
    stack = []  # (depth, domain)
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            rank_code, taxid, name = parts[3], parts[4].strip(), parts[5]
            depth = (len(name) - len(name.lstrip(' '))) // 2
            while stack and stack[-1][0] >= depth:
                stack.pop()
            if rank_code == 'D':
                dom = name.strip()
                current = dom if dom in KEEP_DOMAINS else None
            else:
                current = stack[-1][1] if stack else None
            stack.append((depth, current))
            domain_by_taxid[taxid] = current
    return domain_by_taxid

def parse_bracken(path, domain_by_taxid):
    """flat .bracken file -> ({taxon_name: new_est_reads}, dropped_reads), Bacteria/Archaea only"""
    out = {}
    dropped = 0
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        try:
            name_i = header.index('name')
            taxid_i = header.index('taxonomy_id')
            reads_i = header.index('new_est_reads')
        except ValueError:
            return out, dropped
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= max(name_i, taxid_i, reads_i):
                continue
            try:
                reads = int(parts[reads_i])
            except ValueError:
                continue
            taxid = parts[taxid_i].strip()
            if domain_by_taxid.get(taxid) not in KEEP_DOMAINS:
                dropped += reads
                continue
            name = parts[name_i].strip().replace('"', '').replace("'", '').replace('#', '')
            out[name] = out.get(name, 0) + reads
    return out, dropped

for study, cfg in STUDIES.items():
    print(f"\n{'='*20} {study}")
    kraken_dir = os.path.join(DATA_ROOT, 'work', 'kraken_wgs', study)
    outdir = os.path.join(DATA_ROOT, 'input', cfg['repodir'], 'Kraken2')
    os.makedirs(outdir, exist_ok=True)

    # One domain map per sample (rank-independent — built once from the sample's single
    # .k2report, reused for all 6 ranks' .bracken files below), matching the 16S build
    # step's "drop non-Bacteria/Archaea before normalization" behavior.
    k2reports = sorted(glob.glob(os.path.join(kraken_dir, '*.k2report')))
    domain_maps = {os.path.basename(p)[:-len('.k2report')]: parse_k2report_domains(p) for p in k2reports}
    print(f"k2reports: {len(k2reports)} (expect {cfg['n']})")

    for rank in RANKS:
        files = sorted(glob.glob(os.path.join(kraken_dir, f'*.{rank}.bracken')))
        print(f"{rank}: {len(files)} .bracken files (expect {cfg['n']})")
        if not files:
            print("  none found, skipping (run 03b/04b first)")
            continue

        samples = [os.path.basename(p)[:-len(f'.{rank}.bracken')] for p in files]
        parsed, dropped = {}, 0
        for p, s in zip(files, samples):
            parsed[s], d = parse_bracken(p, domain_maps.get(s, {}))
            dropped += d
        print(f"  non-bacterial/archaeal reads dropped before normalization: {dropped:,}")

        taxa = sorted({t for s in samples for t in parsed[s]})
        totals = {s: sum(parsed[s].values()) for s in samples}
        nonzero = [totals[s] for s in samples if totals[s] > 0]
        mean_total = statistics.mean(nonzero) if nonzero else 0.0

        outfile = os.path.join(outdir, f"{cfg['prefix']}_taxaCount_norm_Log10_{rank}.tsv")
        with open(outfile, 'w') as fh:
            fh.write('SAMPLE_ID\t' + '\t'.join(taxa) + '\n')
            for s in sorted(samples):
                tot = totals[s]
                vals = []
                for t in taxa:
                    c = parsed[s].get(t, 0)
                    v = math.log10(c / tot * mean_total + 1) if tot > 0 and c > 0 else 0.0
                    vals.append(repr(v) if v else '0.0')
                fh.write(s + '\t' + '\t'.join(vals) + '\n')
        print(f"  wrote {outfile}  ({len(taxa)} taxa, mean depth {mean_total:.0f})")

    # metadata: always from this repo's cohorts/<study>/, per the NEEDS_METADATA.md files
    src = os.path.join(COHORTS, study, 'metaData.txt')
    dst_dir = os.path.join(DATA_ROOT, 'input', cfg['repodir'], 'Metadata')
    os.makedirs(dst_dir, exist_ok=True)
    if os.path.exists(src):
        shutil.copy2(src, os.path.join(dst_dir, 'metaData.txt'))
        print(f"  metadata copied from {src}")
        needs = os.path.join(COHORTS, study, 'NEEDS_METADATA.md')
        if os.path.exists(needs):
            print(f"  !! {needs} exists — this cohort's Arm/time/ID columns may still be")
            print(f"     UNKNOWN placeholders. Resolve before trusting any MLM output below.")
    else:
        print(f"  !! no metadata found at {src}")

print("\nDone.")
PYEOF

python3 "$DATA_ROOT/work/build_tables_wgs.py" "$DATA_ROOT" "$COHORTS"
rc=$?

echo
echo "=== regenerated WGS input tree"
find "$DATA_ROOT/input" -path "*_WGS/*" -type f 2>/dev/null | sort | sed 's/^/   /'
exit $rc
