#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 12:00:00
#SBATCH --mem 1G


if [ $# -ne 3 ]; then
    echo "usage: $0 min_depth input.txt output.bed"
    exit 1
fi

mindp=$1
intxt=$2
outbed=$3

if [ -f $outbed ]; then
    echo "output file $outbed already exists, please delete it first"
    exit 1
fi


awk 'BEGIN { OFS="\t"} { if ($1 != "chr" && $3 >= '$mindp') print $1, $2, $2+1; }' $intxt \
    | bedtools merge > $outbed
