#!/bin/bash
#SBATCH -c 1
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=8G

if [ $# -ne 3 ]; then
    echo "usage: $0 sc_dir output_dir n"
    echo "sc_dir must contain training data named 'chrall.rda'"
    exit 1
fi

indir=$1
outdir=$2
n=$3

hsnp_positions.R $indir/training_chrall.rda $outdir/hsnp_positions.txt $n
