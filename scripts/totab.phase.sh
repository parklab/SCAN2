#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: $0 vcf"
    exit 1
fi

vcf=$1

awk 'BEGIN {
    OFS = "\t";
    print "chr", "pos", "ref", "alt", "gt", "hap1", "hap2", "dp", "phgt";
}
{
    if ($1 !~ /^#/) {
        split($10, gtstring, ":");
        gt = gtstring[1];
        if (gt !~ /\./ && $9 ~ ":AD:") {
            dp = gtstring[3];
            split(gtstring[2], adp, ",");
            adp1 = adp[1];
            adp2 = adp[2];
    
            if ($11 == "0|1") {
                print $1, $2, $4, $5, gt, adp1, adp2, dp, $11;
            } else if ($11 == "1|0") {
                print $1, $2, $4, $5, gt, adp2, adp1, dp, $11;
            } else {
                print "ERROR: NOT GOALPOST SITE"
            }
        }
    }
}' $vcf
