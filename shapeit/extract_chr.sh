#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: extract_chr.sh chr"
    exit 1
fi

chr=$1
SHAPEIT_ROOT=/n/data1/hms/dbmi/park/jluquette/pellman/lodato/shapeit


echo "Converting from SHAPEIT format to VCF."
$SHAPEIT_ROOT/bin/shapeit -convert \
    --input-haps shapeit_output_$chr/1000g.cdt.nomulti.chr$chr.phased \
    --output-vcf shapeit_output_$chr/1000g.cdt.nomulti.chr$chr.phased.vcf

echo "Filtering to het sites."
awk '{ if ($10 == "1|0" || $10 == "0|1" || $1 ~ /^#/) print $0; }' \
    shapeit_output_$chr/1000g.cdt.nomulti.chr$chr.phased.vcf > \
    1000g.cdt.nomulti.chr$chr.phased.hets.vcf
