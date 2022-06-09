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

awk '
    function printent(s) {
        split(s, x, ":");
        if (length(x)==1 || x[2] == ".") {
            ad[1]=0;
            ad[2]=0;
        } else {
            split(x[2], ad, ",");
        }
        # Add up contributions from all alleles
        #printf "\n";
        #printf "x[1]=%d, length(ad) = %s\n", x[1], length(ad);
        nalt = 0;
        for (j = 2; j <= length(ad); ++j) {
            nalt = nalt + ad[j];
            #printf "i=%d, ad[j]=%d, nalt=%d\n", j, ad[j], nalt;
        }
        #printf "\n";
        printf "%s\t%s\t%s", x[1], ad[1], nalt;
    }
    {
        split($5, alleles, ",");
        nalleles = length(alleles);
        if ($0 !~ /#/) {
            printf "%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5, nalleles;
            for (i = 10; i <= NF; ++i) {
                printf "\t";
                printent($i);
            }
            printf "\n";
        } else {
            if ($1 ~ /#CHROM/) {
                printf "chr\tpos\tdbsnp\trefnt\taltnt\tnalleles";
                for (i = 10; i <= NF; ++i)
                    printf "\t%s\tref\talt", $i;
                printf "\n";
            }
        }
    }' $INPUT > $OUTPUT
