#!/bin/bash

if [ $# -ne 4 ]; then
    echo "usage: $0 all_regions.bed bad_regions.bed new_window_size out.bed"
    echo "IMPORTANT: bad_regions.bed must be an exact subset of all_regions.bed!"
    echo "If not, some regions in all_regions.bed will be changed and thus trigger reruns."
    echo "This script requires bedtools."
    exit 1
fi

allbed=$1
badbed=$2
window_size=$3
outbed=$4

if [ -f $outbed ]; then
    echo "output file $outbed already exists, please delete it first"
    exit 1
fi

tmpgenome=$(mktemp genome_file.XXXXX.bed)
cut -f1 $allbed | uniq > $tmpgenome

(bedtools subtract -a $allbed -b $badbed ;
    bedtools makewindows -b $badbed -w $window_size ) \
| bedtools sort -g $tmpgenome -i /dev/stdin \
> $outbed

rm $tmpgenome
