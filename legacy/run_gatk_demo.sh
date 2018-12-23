#!/bin/bash
#SBATCH -c 4
#SBATCH -t 12:00:00
#SBATCH -p park
#SBATCH --array=1-396
#SBATCH --mem-per-cpu=5G

if [ $# -lt 5 ]; then
    echo "usage: $0 nthreads mmq outdir bam1 bam2 [bam3 ... ]"
    exit 1
fi

ncores=$1
mmq=$2
outdir=$3
shift 3
bams=$@

mkdir -p $outdir

GATK=$GATK_PATH/gatk.jar
HG19=$GATK_PATH/human_g1k_v37_decoy.fasta
DBSNP=$GATK_PATH/dbsnp.vcf

# If on a SLURM system and submitted as an array, parallelize over the
# regions specified in regions.txt.
vcfid=
region_flag=
if [ ! -z ${SLURM_ARRAY_TASK_ID+x} ]; then
    region_flag="-L $(awk "NR == $SLURM_ARRAY_TASK_ID" regions.txt)"
    echo $region_flag
    vcfid=".${SLURM_ARRAY_TASK_ID}"
fi
region_flag="-L 22"  # Demo only: run on chr22

java -Xms20G -Xmx20G -jar $GATK \
        -nct $ncores \
        -mmq $mmq \
        -T HaplotypeCaller \
        -R $HG19 \
        $(for b in $bams; do echo -I $b; done) \
        --dontUseSoftClippedBases \
        --dbsnp $DBSNP \
        -l INFO \
        -rf BadCigar $region_flag \
        -o $outdir/hc_raw${vcfid}.mmq${mmq}.vcf
