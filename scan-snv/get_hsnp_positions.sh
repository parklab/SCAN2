#!/bin/bash
#SBATCH -c 1
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=8G

if [ $# -ne 2 ]; then
    echo "usage: $0 sc_dir output_dir"
    echo "sc_dir must contain training data named 'chrall.rda'"
    exit 1
fi

indir=$1
outdir=$2

script=/home/ljl11/balance/scan-snv/germline_positions.R

Rscript $script $indir/chrall.rda $outdir/germline_positions.txt
