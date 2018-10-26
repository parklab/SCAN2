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


# Version string is, e.g.:
#   Version: 1.3.1 (using htslib 1.3.1)
# or
#   Version: 0.1.18 (r982:295)
module load gcc
module load samtools/1.3.1
stversion=$(samtools |& grep Version | cut -f2 -d\ )
if [ $stversion != "1.3.1" ]; then
    echo "wrong version of samtools ($stversion), must use 1.3.1"
    exit 1
fi

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
        x=$(samtools view $bam $chrom:$pos-$pos| cut -f6)
        echo "$chrom"$'\t'"$pos"$'\t'$x
    done
} < $infile >> $outfile
