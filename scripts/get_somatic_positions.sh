#!/bin/bash
#SBATCH -c 1
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=12G

if [ $# -ne 5 ]; then
    echo "usage: $0 mmq60.table mmq1.table sc_sample bulk_sample output_dir"
    echo "sc_dir must match the sample name in the mmqXX tables"
    exit 1
fi

mmq60=$1
mmq1=$2
sc=$3
bulk=$4
outdir=$5
ssdir=$outdir

mkdir -p $ssdir

somatic_positions.R $mmq60 $mmq1 $sc $bulk $ssdir/somatic_positions.txt
