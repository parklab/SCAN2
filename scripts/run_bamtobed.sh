#!/bin/bash

#SBATCH -t 12:00:00
#SBATCH -p short
#SBATCH --mem=2G

if [ $# -ne 2 ]; then
    echo "usage: $0 in.bam out.bed"
    exit 1
fi

bam=$1
outbed=$2

if [ -f $outbed ]; then
    echo "output file $outbed already exists, please delete it first"
    exit 1
fi

bedtools bamtobed -i $bam | sed -e 's/^/chr/' > $outbed
