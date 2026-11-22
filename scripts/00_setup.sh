#!/bin/bash
# 00_setup.sh — one-time environment and reference-data setup for the RYGB reproduction.
# Run on a LOGIN or DTN node (it downloads a lot; nothing here needs a compute node):
#     bash 00_setup.sh 2>&1 | tee ~/rygb_setup.log
#
# Layout:  $SW_ROOT (home, 500G, no purge)  = conda envs, singularity images, databases
#          $DATA_ROOT (scratch, 5T)         = raw reads, intermediates, pipeline output
# Every step is idempotent: re-running skips what already exists.

set -uo pipefail

SW_ROOT=${SW_ROOT:-$HOME/rygb}
DATA_ROOT=${DATA_ROOT:-/scratch/$USER/rygb}
ENV_DIR=$SW_ROOT/envs
IMG_DIR=$SW_ROOT/images
DB_DIR=$SW_ROOT/db

# Phase 1 = 16S only. Set to 1 when you move on to the metagenomic half.
SETUP_METAGENOMICS=${SETUP_METAGENOMICS:-0}

# Kraken2 index. The authors' custom DB (/nobackup/afodor_research/...) died with Copperhead.
# Closest public stand-in by date; check https://benlangmead.github.io/aws-indexes/k2 if this 404s.
K2_URL=${K2_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20201202.tar.gz}

sec(){ printf '\n\033[1m==== %s\033[0m\n' "$1"; }
have(){ [ -e "$1" ]; }

sec "0. Directories"
mkdir -p "$ENV_DIR" "$IMG_DIR" "$DB_DIR" \
         "$DATA_ROOT"/{raw,work,input,pipeline,logs}
mkdir -p "$DATA_ROOT"/raw/{BS16S,Assal,Afshar,Ilhan,BSmeta,Palleja}
echo "software: $SW_ROOT"
echo "data:     $DATA_ROOT"
df -hP "$HOME" "$DATA_ROOT" | sed 's/^/  /'

sec "1. Modules"
module load mambaforge/23.11
module load singularity/4.4.1
module load sra-tools/3.1.0
module load kraken2/2.1.3
module list 2>&1 | sed 's/^/  /'

export CONDA_PKGS_DIRS=$SW_ROOT/conda_pkgs
export SINGULARITY_CACHEDIR=$SW_ROOT/singularity_cache
mkdir -p "$CONDA_PKGS_DIRS" "$SINGULARITY_CACHEDIR"

CH="-c conda-forge -c bioconda -c defaults --override-channels"

sec "2. Conda env: rygb-16s  (DADA2 + cutadapt)"
if have "$ENV_DIR/rygb-16s/bin/R"; then
  echo "  exists, skipping"
else
  mamba create -y -p "$ENV_DIR/rygb-16s" $CH \
    r-base=4.0.3 bioconductor-dada2=1.16.0 bioconductor-shortread \
    cutadapt=3.4 pigz
fi

if [ "$SETUP_METAGENOMICS" = "1" ]; then
  sec "3a. Conda env: rygb-kneaddata (QC + human-read removal, plus Bracken for the taxonomy-only WGS comparison track — see section 6b below)"
  if have "$ENV_DIR/rygb-kneaddata/bin/kneaddata"; then
    echo "  exists, skipping"
  else
    mamba create -y -p "$ENV_DIR/rygb-kneaddata" $CH \
      kneaddata=0.12.0 bowtie2 trimmomatic bracken=2.9 || \
      echo "  !! solve failed — tell me the error; bracken=2.9 may conflict with kneaddata=0.12.0's pinned deps"
  fi
  if have "$ENV_DIR/rygb-kneaddata/bin/bracken"; then
    echo "  bracken present"
  else
    echo "  !! bracken missing from an existing rygb-kneaddata env — installing into it"
    mamba install -y -p "$ENV_DIR/rygb-kneaddata" $CH bracken=2.9
  fi

  # kneaddata's --trimmomatic expects a DIRECTORY it globs "trimmomatic*" in via Python's
  # os.listdir() + fnmatch, taking found_paths[0] — bioconda's trimmomatic package ships
  # both a Python wrapper script named "trimmomatic" and the real "trimmomatic.jar" in the
  # same share/ dir, and os.listdir() order (not alphabetical) can return the wrapper
  # first, so kneaddata tries to `java -jar` the wrapper script itself and fails ("Invalid
  # or corrupt jarfile"). Confirmed via kneaddata's own source
  # (kneaddata/utilities.py:find_dependency). Fix: a clean directory containing ONLY the
  # real jar (symlinked), so there's no ambiguity — 03b_kneaddata.sbatch points
  # --trimmomatic at this instead of the raw share/ dir.
  TRIMMOMATIC_SHARE=$(ls -d "$ENV_DIR/rygb-kneaddata"/share/trimmomatic*/ 2>/dev/null | head -1)
  TRIMMOMATIC_JAR_DIR=$DB_DIR/trimmomatic_jar_only
  if [ -n "$TRIMMOMATIC_SHARE" ] && [ -e "$TRIMMOMATIC_SHARE/trimmomatic.jar" ]; then
    if have "$TRIMMOMATIC_JAR_DIR/trimmomatic.jar"; then
      echo "  trimmomatic_jar_only/ already set up"
    else
      mkdir -p "$TRIMMOMATIC_JAR_DIR"
      ln -sf "$TRIMMOMATIC_SHARE/trimmomatic.jar" "$TRIMMOMATIC_JAR_DIR/trimmomatic.jar"
      echo "  built $TRIMMOMATIC_JAR_DIR (symlink to $TRIMMOMATIC_SHARE/trimmomatic.jar)"
    fi
  else
    echo "  !! no trimmomatic.jar found under $ENV_DIR/rygb-kneaddata/share — check the env"
  fi

  # kneaddata's bundled kneaddata/trf_parallel.py has a real upstream bug: run_trf()
  # always splits its input into exactly --nproc chunks via integer division
  # (lines_per_file = int(total_lines/nproc)) regardless of actual sequence count. For a
  # small input — routine for kneaddata's own "trimmed.single.*" orphaned-read files —
  # most chunks end up empty, and the original code runs trf on every one of them
  # unconditionally, so trf fails on each empty chunk ("Could not load sequence. Empty
  # file or bad format."), which kneaddata then reports as "Subprocess reported error"
  # and the whole sample fails. Confirmed via direct reproduction (a 1-sequence input and,
  # separately, a 50-sequence input at nproc=16, both left empty trailing chunks) before
  # this was ever patched live in a running env — codified here so a fresh/rebuilt env
  # gets the fix automatically instead of silently losing it (this previously lived only
  # as an unpatched live edit under $ENV_DIR, which is gitignored and was never tracked).
  # Glob rather than hardcode the python3.12 directory name: that version is whatever
  # mamba's solver happened to pick for kneaddata=0.12.0's dependency chain this time,
  # not something pinned explicitly above — a re-run on a different day/channel state,
  # or a fresh machine, is not guaranteed to land on the same python minor version this
  # session actually got, and a hardcoded path would then silently miss the installed
  # file entirely (falling into the "not found" branch below with a misleading "will be
  # patched next run" message, since the path would stay wrong forever, not just absent).
  TRF_PARALLEL=$(ls "$ENV_DIR"/rygb-kneaddata/lib/python3.*/site-packages/kneaddata/trf_parallel.py 2>/dev/null | head -1)
  if [ -n "$TRF_PARALLEL" ] && [ -e "$TRF_PARALLEL" ]; then
    if grep -q "RYGB patch" "$TRF_PARALLEL" 2>/dev/null; then
      echo "  trf_parallel.py already patched"
    else
      cat > "$TRF_PARALLEL" <<'TRFPY'
#!/usr/bin/env python

"""
KneadData parallel TRF

This script wraps bowtie2 to allow for discordant paired alignments output to fastq.

Dependencies: TRF

Copyright (c) 2017 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import sys
import os

import argparse

import string
import tempfile
import subprocess

try:
    from kneaddata import utilities
except ImportError:
    sys.exit("Please install kneaddata")

def parse_arguments(args):
    """
    Parse the arguments from the user
    """

    parser = argparse.ArgumentParser(
    description= "Kneaddata trf parallel\n",
    formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--input",
        help="the fasta file of reads",
        required=True)
    parser.add_argument(
        "--output",
        help="the trf output file",
        required=True)
    parser.add_argument(
        "--trf-path",
        help="full path to the trf executable",
        required=True)
    parser.add_argument(
        "--trf-options",
        help="space delimited list of trf options",
        required=True)
    parser.add_argument(
        "--nproc",
        default=1,
        type=int,
        help="total number of processes to use")

    return parser.parse_args()


def run_trf(input,trf_path,trf_options,nproc,output,verbose=True):
    """ Run trf with the options provided """

    tempfile_list=[]
    datfile_list=[]
    commands=[]

    # check for one process and if so just run trf directly
    if nproc == 1:
        commands.append([[trf_path,input]+trf_options.split(" "),"trf",[input],[output],output])

        utilities.start_processes(commands,nproc,verbose)
    else:
        # get the total number of reads
        total_lines=0
        with open(input) as file_handle:
            for line in file_handle:
                total_lines+=1

        # split the input into multiple files and run in parallel
        for i in range(int(nproc)):
            file_out, new_file = tempfile.mkstemp(prefix=os.path.basename(output)+'_'+str(i)+'_temp_trf_output',dir=os.path.dirname(output))
            os.close(file_out)
            tempfile_list.append(new_file)
            datfile_list.append(new_file+".".join(trf_options.split("-")[0].split(" "))+"dat")

        # write the input file into all temp output files
        output_file_number=0
        lines_per_file = int(total_lines/int(nproc))
        lines_written=0
        file_handle_write=None
        for read_line in utilities.read_file_n_lines(input,2):
            if not file_handle_write:
                file_handle_write = open(tempfile_list[output_file_number],"wt")
            file_handle_write.write("".join(read_line))

            lines_written+=2
            if lines_written > lines_per_file:
                file_handle_write.close()
                lines_written=0
                output_file_number+=1
                file_handle_write = open(tempfile_list[output_file_number],"wt")

        file_handle_write.close()

        # RYGB patch: when total_lines doesn't divide evenly across nproc chunks (routine
        # whenever the input is small — e.g. kneaddata's "trimmed.single.*" orphaned-read
        # files), the write loop above leaves one or more trailing tempfile_list entries
        # completely empty (0 bytes, never opened by the loop, or opened-but-never-written
        # as the final "next" file). Running trf on an empty file always fails ("Could not
        # load sequence. Empty file or bad format."), and the original code ran it
        # unconditionally on every one of the nproc chunks regardless of content. Confirmed
        # via direct reproduction (both a 1-sequence input and, separately, a 50-sequence
        # input with nproc=16, both left empty trailing chunks) before patching. Fix: only
        # run trf on chunks that actually received content — skip the rest, they have
        # nothing to classify anyway.
        nonempty = [(i, t, d) for i, (t, d) in enumerate(zip(tempfile_list, datfile_list))
                    if os.path.getsize(t) > 0]

        # run commands
        for i, temp_in, temp_out in nonempty:
            trf_command=[trf_path,temp_in]+trf_options.split(" ")
            commands.append([trf_command,"trf{}".format(i),[temp_in],[temp_out],temp_out])

        utilities.start_processes(commands,nproc,verbose)

        # datfile_list only needs to cover chunks actually run, for the merge step below
        datfile_list = [d for _, _, d in nonempty]

        # merge all of the outputs to the final output file
        with open(output,"w") as file_write:
            for datfile in datfile_list:
                with open(datfile) as file_read:
                    for line in file_read:
                        file_write.write(line)

        # remove temp files
        for filename in tempfile_list+datfile_list:
            try:
                os.remove(filename)
            except EnvironmentError:
                print("Unable to remove temp file: " + filename)

def main():
    # parse the command line arguments
    args = parse_arguments(sys.argv)

    # remove extra quotes if included
    args.trf_options=args.trf_options.replace("'","")

    # add exe name if not included
    if not args.trf_path.endswith("trf"):
        if not args.trf_path.endswith("/"):
            args.trf_path+="/"
        args.trf_path+="trf"

    run_trf(args.input,args.trf_path,args.trf_options,args.nproc,args.output)

if __name__ == "__main__":
    main()
TRFPY
      rm -f "$(dirname "$TRF_PARALLEL")/__pycache__/trf_parallel.cpython"*.pyc 2>/dev/null
      echo "  patched $TRF_PARALLEL (empty-chunk trf fix)"
    fi
  else
    echo "  !! kneaddata/trf_parallel.py not found under $ENV_DIR/rygb-kneaddata — kneaddata env not built yet? (will be patched next run)"
  fi

  sec "3b. Conda env: rygb-humann2  (Python 2.7 — the fragile one)"
  # Only needed for pathway/functional profiling, which is OUT OF SCOPE for the current
  # taxonomy-only WGS comparison (see section 6b below) — this block and the HUMAnN2
  # database download in section 7 are left as-is from the original setup but are not
  # required by scripts/03b_kneaddata.sbatch, 04b_kraken2_wgs.sbatch, or 05b_build_tables_wgs.sbatch.
  # HUMAnN2 2.8.1 needs diamond 0.8.36 exactly; the cluster's diamond module (2.x) will not work.
  if have "$ENV_DIR/rygb-humann2/bin/humann2"; then
    echo "  exists, skipping"
  else
    mamba create -y -p "$ENV_DIR/rygb-humann2" $CH \
      python=2.7 humann2=2.8.1 metaphlan2=2.7.7 diamond=0.8.36 bowtie2=2.3.5.1 || \
      echo "  !! solve failed — tell me the error; fallback is a biobakery container"
  fi
fi

sec "4. Singularity images (the authors' R environments)"
for img in rbase nlme r-vegan comparestudies complexheatmap pathogens siamcat fig-perform; do
  if have "$IMG_DIR/$img.sif"; then
    echo "  $img.sif exists"
  else
    echo "  pulling asorgen/$img:v1"
    singularity pull "$IMG_DIR/$img.sif" "docker://asorgen/$img:v1" \
      || echo "  !! pull failed for $img — note which ones fail"
  fi
done

sec "5. SILVA 132 (DADA2 training sets)"
mkdir -p "$DB_DIR/silva132"
for f in silva_nr_v132_train_set.fa.gz silva_species_assignment_v132.fa.gz; do
  if have "$DB_DIR/silva132/$f"; then
    echo "  $f exists"
  else
    curl -fL --retry 3 -o "$DB_DIR/silva132/$f" \
      "https://zenodo.org/record/1172783/files/$f?download=1" \
      || echo "  !! download failed: $f"
  fi
done
ls -lh "$DB_DIR/silva132" | sed 's/^/  /'

sec "6. Kraken2 index"
if have "$DB_DIR/kraken2/hash.k2d"; then
  echo "  exists, skipping"
  du -sh "$DB_DIR/kraken2" | sed 's/^/  /'
else
  echo "  source: $K2_URL"
  code=$(curl -s -o /dev/null -m 30 -w '%{http_code}' -I "$K2_URL")
  if [ "$code" != "200" ]; then
    echo "  !! HTTP $code — index not there. Pick another from"
    echo "     https://benlangmead.github.io/aws-indexes/k2 and rerun with K2_URL=..."
  else
    mkdir -p "$DB_DIR/kraken2"
    curl -fL --retry 3 "$K2_URL" | tar -xz -C "$DB_DIR/kraken2"
    du -sh "$DB_DIR/kraken2" | sed 's/^/  /'
    ls -lh "$DB_DIR/kraken2"/*.k2d | sed 's/^/  /'   # hash.k2d size sets the job's --mem
  fi
fi

if [ "$SETUP_METAGENOMICS" = "1" ]; then
  sec "6b. Bracken build (companion to the Kraken2 index, WGS taxonomy track only)"
  # bracken-build re-uses the Kraken2 index built in section 6 — it does NOT need its own
  # database download, just a one-time k-mer distribution pass over that same index, once
  # per read length actually used by 04b_kraken2_wgs.sbatch's per-study READLEN map: 150bp
  # (Penney, Oseberg, COCKTAIL_GBP — all ~150bp/NovaSeq 2x150) and 300bp (Hayashi_LSG,
  # MiSeq 2x300 — far enough from 150bp to need its own distribution; Bracken's own
  # guidance is ~50bp tolerance).
  BRACKEN_READLENS=${BRACKEN_READLENS:-"150 300"}
  for BRACKEN_READLEN in $BRACKEN_READLENS; do
    if have "$DB_DIR/kraken2/database${BRACKEN_READLEN}mers.kmer_distrib"; then
      echo "  exists for ${BRACKEN_READLEN}bp, skipping"
    else
      "$ENV_DIR/rygb-kneaddata/bin/bracken-build" -d "$DB_DIR/kraken2" \
        -t "$(nproc)" -l "$BRACKEN_READLEN" -k 35 \
        || echo "  !! bracken-build failed for ${BRACKEN_READLEN}bp — needs rygb-kneaddata env built first (section 3a)"
    fi
  done

  sec "7. KneadData human index + HUMAnN2 databases (~30 GB)"
  source "$(dirname "$(dirname "$(which mamba)")")/etc/profile.d/conda.sh"
  if have "$ENV_DIR/rygb-kneaddata/bin/kneaddata_database"; then
    mkdir -p "$DB_DIR/kneaddata_human"
    conda run -p "$ENV_DIR/rygb-kneaddata" \
      kneaddata_database --download human_genome bowtie2 "$DB_DIR/kneaddata_human"
  fi
  if have "$ENV_DIR/rygb-humann2/bin/humann2_databases"; then
    mkdir -p "$DB_DIR/humann2"
    conda run -p "$ENV_DIR/rygb-humann2" humann2_databases --download chocophlan full "$DB_DIR/humann2"
    conda run -p "$ENV_DIR/rygb-humann2" humann2_databases --download uniref uniref90_diamond "$DB_DIR/humann2"
  fi
fi

sec "8. The repository (analysis scripts + metadata we reuse)"
if have "$DATA_ROOT/RYGB_IntegratedAnalysis2020/README.md"; then
  echo "  already cloned"
else
  git clone https://github.com/FarnazFouladi/RYGB_IntegratedAnalysis2020.git \
    "$DATA_ROOT/RYGB_IntegratedAnalysis2020" || echo "  !! clone failed"
fi

sec "9. Provenance"
{
  date
  echo "SW_ROOT=$SW_ROOT  DATA_ROOT=$DATA_ROOT"
  echo "kraken2: $(kraken2 --version 2>&1 | head -1)"
  echo "sra-tools: $(fasterq-dump --version 2>&1 | tr -d '\n')"
  echo "singularity: $(singularity --version)"
  echo "K2_URL=$K2_URL"
  for e in "$ENV_DIR"/*/; do echo "env: $e"; done
  for i in "$IMG_DIR"/*.sif; do echo "img: $i"; done
} > "$SW_ROOT/versions.txt"
cat "$SW_ROOT/versions.txt" | sed 's/^/  /'

sec "Done — report back:"
echo "  * which singularity pulls failed, if any"
echo "  * the hash.k2d size from section 6"
echo "  * any conda solve errors"
