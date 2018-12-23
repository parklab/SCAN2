#!/bin/bash
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -p park
#SBATCH --array=1-396
#SBATCH --mem-per-cpu=10G

if [ $# -lt 4 ]; then
    echo "usage: $0 mmq outdir bam1 bam2 [bam3 ... ]"
    echo "NOTE: if submitting to SLURM, nthreads should match the SBATCH"
    echo "header in this file."
    exit 1
fi

mmq=$1
outdir=$2
shift 2
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

java -Xms8G -Xmx8G -jar $GATK \
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
