#!/bin/bash

if [ $# -ne 2 ]; then
    echo "usage: $0 input_depth_of_coverage_table output"
    exit 1
fi

INPUT=$1
OUTPUT=$2

if [ -f $OUTPUT ] && [ $OUTPUT != '/dev/stdout' ]; then
    echo "output file $OUTPUT already exists, please delete it first"
    exit 1
fi

awk 'BEGIN { OFS="\t"; } {
    if ($0 !~ /^Locus/) {
        split($1, firstcol, ":")
        printf "%s\t%s", firstcol[1], firstcol[2];
            for (i = 4; i <= NF; ++i) {
            printf "\t%d", $i;
        }
        printf "\n";
    } else {
        printf "#chr\tpos";
        for (i = 4; i <= NF; ++i) {
            sub(/^Depth_for_/, "", $i);
            printf "\t%s", $i;
        }
        printf "\n";
    }
}' $INPUT > $OUTPUT
