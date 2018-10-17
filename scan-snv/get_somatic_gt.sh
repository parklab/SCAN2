#!/bin/bash
#SBATCH -c 1
#SBATCH -t 6:00:00
#SBATCH -A park_contrib
#SBATCH -p park
#SBATCH --mem-per-cpu=12G

if [ $# -ne 6 ]; then
    echo "usage: $0 mmq60.tab mmq1.tab sc.sample bulk.sample scan_snv_dir"
    exit 1
fi

module load gcc
module load R/3.3.3

mmq60=$1
mmq1=$2
sc_sample=$3
bulk_sample=$4
dir=$5
fdr=$6

script=/home/ljl11/balance/scan-snv/genotype.R

if [ -f $dir/somatic_gt.rda ]; then
    echo "output file $outfile already exists, please delete it first"
    exit 1
fi

Rscript $script $mmq60 $mmq1 $sc_sample $bulk_sample $dir/somatic_ab.rda $dir/somatic_cigars.tab $dir/hsnp_cigars.tab $dir/somatic_gt.rda $fdr
