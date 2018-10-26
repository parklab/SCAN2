#!/bin/bash

if [ $# -ne 2 ]; then
    echo "usage: extract_chr.sh output_dir chr"
    echo "output_dir should be the same path used for output_dir"
    echo "when running phase_chr.sh."
    exit 1
fi

output_dir=$1
chr=$2
fromsample=`cat samples_to_phase.txt`
SHAPEIT_ROOT=/n/data1/hms/dbmi/park/jluquette/pellman/lodato/shapeit


chrdir=$output_dir/shapeit_output_$chr

echo "Converting from SHAPEIT format to VCF."
$SHAPEIT_ROOT/bin/shapeit -convert \
    --input-haps $chrdir/1000g.cdt.nomulti.chr$chr.phased \
    --output-vcf $chrdir/1000g.cdt.nomulti.chr$chr.phased.vcf

echo "Filtering to het sites."
awk '{ if ($10 == "1|0" || $10 == "0|1" || $1 ~ /^#/) print $0; }' \
    $chrdir/1000g.cdt.nomulti.chr$chr.phased.vcf \
    | sed -e"s/$fromsample/phasedgt/g" > \
      $output_dir/phased_hsnps.chr$chr.vcf
