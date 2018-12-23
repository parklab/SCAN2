#!/bin/bash
#SBATCH -c 1
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 12:00:00
#SBATCH --array=1-22
#SBATCH --mem-per-cpu=8G

if [ $# -lt 3 ]; then
    echo "usage: ./phase_chr.sh vcf output_dir bulk_sample [chromosome]"
    exit 1
fi

vcf=$(realpath $1)
output_dir=$2
bulksample=$3

cd $output_dir

if [ ! -z ${SLURM_ARRAY_TASK_ID+x} ]; then
    chr=$SLURM_ARRAY_TASK_ID
else
    if [ $# != 4 ]; then
        echo "Not in SLURM array; must specify chromosome"
        exit 1
    fi
    chr=$4
fi

# SHAPEIT runs in 2 steps.  The first step prunes out SNPs in the
# reference panel that are problematic, step 2 performs the phasing.
OUTPUT_ROOT=shapeit_output_$chr

mkdir -p $OUTPUT_ROOT

echo "Phasing chromosome=$chr from VCF=$vcf"

tmpfile=`mktemp`

cleanvcf=$OUTPUT_ROOT/in.vcf
java -jar $GATK_PATH/gatk.jar \
    -R $GATK_PATH/human_g1k_v37_decoy.fasta \
    -T SelectVariants \
    -V $vcf \
    -selectType SNP -restrictAllelesTo BIALLELIC -env -trimAlternates \
    -select 'vc.getGenotype("'$bulksample'").isCalled()' \
    -L $chr \
    -o $cleanvcf

echo $bulksample > $tmpfile

if [ $chr != "X" ]; then
    gmap=$REFPANEL_ROOT/genetic_map_chr${chr}_combined_b37.txt
    hap=$REFPANEL_ROOT/1000GP_Phase3_chr${chr}.hap.gz
    leg=$REFPANEL_ROOT/1000GP_Phase3_chr${chr}.legend.gz
else
    gmap=$REFPANEL_ROOT/genetic_map_chr${chr}_nonPAR_combined_b37.txt
    hap=$REFPANEL_ROOT/1000GP_Phase3_chr${chr}_NONPAR.hap.gz
    leg=$REFPANEL_ROOT/1000GP_Phase3_chr${chr}_NONPAR.legend.gz
fi

$SHAPEIT_ROOT/bin/shapeit \
    -check \
    --input-vcf=$cleanvcf \
    -M $gmap \
    --input-ref $hap $leg $REFPANEL_ROOT/1000GP_Phase3.sample \
    --include-ind $tmpfile \
    --output-log $OUTPUT_ROOT/shapeit_check.log

xflag=
if [ $chr == "X" ]; then
    xflag="--chrX"
fi

$SHAPEIT_ROOT/bin/shapeit \
    --input-vcf=$cleanvcf \
    -M $gmap \
    --input-ref $hap $leg $REFPANEL_ROOT/1000GP_Phase3.sample \
    --include-ind $tmpfile \
    --exclude-snp $OUTPUT_ROOT/shapeit_check.snp.strand.exclude \
    $xflag \
    -O $OUTPUT_ROOT/1000g.cdt.nomulti.chr$chr.phased



echo "Converting from SHAPEIT format to VCF."
$SHAPEIT_ROOT/bin/shapeit -convert \
    --input-haps $OUTPUT_ROOT/1000g.cdt.nomulti.chr$chr.phased \
    --output-vcf $OUTPUT_ROOT/1000g.cdt.nomulti.chr$chr.phased.vcf

echo "Filtering to het sites."
awk '{ if ($10 == "1|0" || $10 == "0|1" || $1 ~ /^#/) print $0; }' \
    $OUTPUT_ROOT/1000g.cdt.nomulti.chr$chr.phased.vcf \
    | sed -e"s/$bulksample/phasedgt/g" > \
      phased_hsnps.chr$chr.vcf
