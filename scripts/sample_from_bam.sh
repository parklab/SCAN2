#!/bin/bash

if [ $# -lt 1 ]; then
    echo "usage: $0 bam1 [ bam2 [ bam3 ... [ bamN ] ] ]"
    echo "N.B.: assumes only one sample per BAM file"
    exit 1
fi

bams=$@

for b in $bams; do
    # Extract the sample name, assuming only 1 sample per bam
    sn=$(samtools view -H $b|grep -m1 -oP 'SM:[^ \t]*'|cut -f2 -d\:)
    echo "    --bam $sn $b \\"
    echo "    --sc-sample $sn \\"
done 
