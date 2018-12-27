#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu=8G
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 12:00:00

if [ $# -lt 2 ]; then
    echo "usage: $0 output_vcf vcf1 [ vcf2 ... vcfN ]"
    exit 1
fi

outvcf=$1
shift 1

java -Xmx7G -Xms7G \
    -cp $GATK_PATH/gatk.jar org.broadinstitute.gatk.tools.CatVariants \
    -R $GATK_PATH/human_g1k_v37_decoy.fasta \
    $(for vcf in $@; do echo -V $vcf; done) \
    -out $outvcf \
    -assumeSorted
