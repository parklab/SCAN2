#!/bin/bash
#SBATCH -c 4
#SBATCH -t 12:00:00
#SBATCH -p park
#SBATCH --array=1-396
#SBATCH --mem-per-cpu=5G

if [ $# -lt 4 ]; then
    echo "usage: $0 mmq outdir bam1 bam2 [bam3 ... ]"
    exit 1
fi

mmq=$1
outdir=$2
shift
shift
bams=$@

mkdir -p $outdir

# ------ EDIT THESE VARIABLES ------
# If multiple cores are available, specify the number of cores here.
ncores=8
mem=22G     # If using ncores > 1, increase ~linearly up to ~24G

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

java -Xms${mem} -Xmx${mem} -jar $GATK \
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
