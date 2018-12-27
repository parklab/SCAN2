#!/bin/bash
#SBATCH -c 1
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH --mem-per-cpu=8G
#SBATCH -t 12:00:00

if [ $# -ne 4 ]; then
    echo "usage: $0 single_cell_sample joint_vcf phased_vcf output_dir"
    exit 1
fi

sample=$1
joint_vcf=$2
phased_vcf=$3
output_dir=$4

output_vcf=$output_dir/${sample}_hsnps.vcf
output_tab=$output_dir/${sample}_hsnps.tab
tmp_vcf=$output_dir/tmp.vcf

GATK=$GATK_PATH/gatk.jar
REF=$GATK_PATH/human_g1k_v37_decoy.fasta

java -Xmx6G -Xms6G \
    -jar $GATK -R $REF \
    -T CombineVariants \
    -V $joint_vcf \
    -V $phased_vcf \
    -o $tmp_vcf

java -Xmx6G -Xms6G \
    -jar $GATK -R $REF \
    -T SelectVariants \
    -V $tmp_vcf \
    -sn $sample \
    -sn phasedgt \
    -env -trimAlternates \
    -select 'vc.getGenotype("'$sample'").isCalled()' \
    -select 'vc.getGenotype("phasedgt").isCalled()' \
    -select 'vc.isBiallelic()' \
    -selectType SNP \
    -o $output_vcf

rm $tmp_vcf
rm ${tmp_vcf}.idx

totab.phase.sh $output_vcf $output_tab

torda.R $output_tab $output_dir/training_chr%d
