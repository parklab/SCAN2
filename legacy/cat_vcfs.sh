#!/bin/bash

if [ $# -lt 2 ]; then
    echo "usage: $0 output_vcf vcf1 [ vcf2 ... vcfN ]"
    exit 1
fi

outvcf=$1
shift 1

java -cp $GATK_PATH/gatk.jar org.broadinstitute.gatk.tools.CatVariants \
    -R $GATK_PATH/human_g1k_v37_decoy.fasta \
    $(for vcf in $@; do echo -V $vcf; done) \
    -out $outvcf \
    -assumeSorted
