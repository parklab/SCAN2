#!/bin/bash
#SBATCH -c 1
#SBATCH -t 6:00:00
#SBATCH -A park_contrib
#SBATCH -p park

if [ $# -ne 3 ]; then
    echo "usage: $0 position_file bam outfile"
    exit 1
fi

infile=$1
bam=$2
outfile=$3


if [ -f $outfile ]; then
    echo "output file $outfile already exists, please delete it first"
    exit 1
fi

echo "chrom"$'\t'"pos"$'\t'"cigars" >> $outfile

{
    read # skip header line
    while read line; do
        array=($line)
        chrom=${array[0]}
        pos=${array[1]}
        x=$(samtools view -q 60 $bam $chrom:$pos-$pos| cut -f6)
        echo "$chrom"$'\t'"$pos"$'\t'$x
    done
} < $infile >> $outfile
