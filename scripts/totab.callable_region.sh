#!/bin/bash

if [ $# -ne 2 ]; then
    echo "usage: $0 inputvcf output"
    exit 1
fi

INPUT=$1
OUTPUT=$2

if [ -f $OUTPUT ] && [ $OUTPUT != '/dev/stdout' ]; then
    echo "output file $OUTPUT already exists, please delete it first"
    exit 1
fi

awk 'BEGIN { OFS="\t"; } {
    if ($0 !~ /#/) {
        split($1, firstcol, ":")
        printf "%s\t%s\t%s\t%s\t%s\t%s", firstcol[1], firstcol[2];
            for (i = 4; i <= NF; ++i) {
            print "\t" $i;
        }
        print "\n";
    } else {
        if ($1 ~ /^Locus/) {
            print "#chr\tpos";
            for (i = 4; i <= NF; ++i) {
                print "\t" sub(/^Depth_for_/, "", $i);
            }
            print "\n";
        }
    }
}' $INPUT > $OUTPUT
