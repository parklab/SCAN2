#!/bin/bash
#SBATCH -c 1
#SBATCH -t 6:00:00
#SBATCH -A park_contrib
#SBATCH -p park
#SBATCH --mem-per-cpu=8G

if [ $# -ne 2 ]; then
    echo "usage: $0 dir output_dir"
    exit 1
fi

dir=$1
outdir=$2

infile=$outdir/somatic_positions.txt
training=$dir/training_chrall.rda
fit=$dir/fits.rda
outfile=$outdir/somatic_ab.rda

if [ -f $outfile ]; then
    echo "output file $outfile already exists, please delete it first"
    exit 1
fi

somatic_ab.R $infile $training $fit $outfile
