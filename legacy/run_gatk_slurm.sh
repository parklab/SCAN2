#!/bin/bash
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH --array=1-410
#SBATCH --mem-per-cpu=10G

if [ $# -lt 3 ]; then
    echo "usage: $0 outdir bam1 bam2 [bam3 ... ]"
    exit 1
fi

outdir=$1
shift 1
bams=$@

mkdir -p $outdir

GATK=$GATK_PATH/gatk.jar
HG19=$GATK_PATH/human_g1k_v37_decoy.fasta
DBSNP=$GATK_PATH/dbsnp.vcf
REGIONS=$GATK_PATH/regions.txt

# If on a SLURM system and submitted as an array, parallelize over the
# regions specified in regions.txt.
vcfid=
region_flag=
if [ ! -z ${SLURM_ARRAY_TASK_ID+x} ]; then
    # Array values 1-205 use MMQ=60; 206-410 use MMQ=1.
    # There are only 200 regions in the standard regions file.
    mmq=60
    regionidx=$SLURM_ARRAY_TASK_ID
    if [ $SLURM_ARRAY_TASK_ID -gt 205 ]; then
        mmq=1
        regionidx=$(($regionidx - 205))
    fi

    region_flag="-L $(awk "NR == $regionidx" $REGIONS)"
    echo $region_flag
    vcfid=".$regionidx"
fi


java -Xms8G -Xmx8G -jar $GATK \
        -mmq $mmq \
        -T HaplotypeCaller \
        -R $HG19 \
        $(for b in $bams; do echo -I $b; done) \
        --dontUseSoftClippedBases \
        --dbsnp $DBSNP \
        -l INFO \
        -rf BadCigar $region_flag \
        -o $outdir/hc_raw${vcfid}.mmq${mmq}.vcf
