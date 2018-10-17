#!/bin/bash
#SBATCH -c 1
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 12:00:00
#SBATCH --array=1-22
#SBATCH --mem-per-cpu=8G

if [ $# != 1 ]; then
    echo "usage: ./phase_chr.sh vcf"
    exit 1
fi

vcf=$1
chr=$SLURM_ARRAY_TASK_ID

echo "Phasing chromosome=$chr from VCF=$vcf"

chrvcf=${vcf/.vcf/.chr$chr.vcf}

# Split by chromosome.  I think SHAPEIT can do this via some command
# line arguments, but was too lazy to find them.
awk "{ if (\$1 == \"$chr\" || \$0 ~ /^#/) print \$0; }" $vcf > $chrvcf


# SHAPEIT runs in 2 steps.  The first step prunes out SNPs in the
# reference panel that are problematic, step 2 performs the phasing.
SHAPEIT_ROOT=/n/data1/hms/dbmi/park/jluquette/pellman/lodato/shapeit
REFPANEL_ROOT=$SHAPEIT_ROOT/ALL.integrated_phase1_SHAPEIT_16-06-14.nosing
OUTPUT_ROOT=./shapeit_output_$chr

mkdir -p $OUTPUT_ROOT

$SHAPEIT_ROOT/bin/shapeit \
    -check \
    --input-vcf=$chrvcf \
    -M $REFPANEL_ROOT/genetic_map_chr${chr}_combined_b37.txt \
    --input-ref $REFPANEL_ROOT/ALL.chr$chr.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.haplotypes.gz \
        $REFPANEL_ROOT/ALL.chr$chr.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz \
        $REFPANEL_ROOT/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample \
    --include-ind samples_to_phase.txt \
    --output-log $OUTPUT_ROOT/shapeit_check.log

$SHAPEIT_ROOT/bin/shapeit \
    --input-vcf=$chrvcf \
    -M $REFPANEL_ROOT/genetic_map_chr${chr}_combined_b37.txt \
    --input-ref $REFPANEL_ROOT/ALL.chr$chr.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.haplotypes.gz \
        $REFPANEL_ROOT/ALL.chr$chr.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz \
        $REFPANEL_ROOT/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample \
    --include-ind samples_to_phase.txt \
    --exclude-snp $OUTPUT_ROOT/shapeit_check.snp.strand.exclude \
    -O $OUTPUT_ROOT/1000g.cdt.nomulti.chr$chr.phased
