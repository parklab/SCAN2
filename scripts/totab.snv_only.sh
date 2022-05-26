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
        printf "%s\t%s\t%s", x[1], ad[1], ad[2];
    }
    {
        if ($0 !~ /#/ && length($4) == 1) { # && length($5) == 1) {
            n = split($8, info, ";");
            mq = "NA";
            mqrs = "NA";
            for (i = 0; i < n; ++i) {
                if (info[i] ~ /^MQ=/) {
                    split(info[i], tmp, "=");
                    mq = tmp[2];
                }
                if (info[i] ~ /^MQRankSum=/) {
                    split(info[i], tmp, "=");
                    mqrs = tmp[2];
                }
            }
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5, mq, mqrs;
            for (i = 10; i <= NF; ++i) {
                printf "\t";
                printent($i);
            } printf "\n";
        } else {
            if ($1 ~ /#CHROM/) {
                printf "chr\tpos\tdbsnp\trefnt\taltnt\tmq\tmqrs";
                for (i = 10; i <= NF; ++i)
                    printf "\t%s\tref\talt", $i;
                printf "\n";
            }
        }
    }' $INPUT > $OUTPUT
