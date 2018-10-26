#!/bin/bash
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -p short
#SBATCH --array=1-396
#SBATCH --mem-per-cpu=10G

if [ $# -ne 1 ]; then
    echo "usage: $0 mmq"
    exit 1
fi

mmq=$1
slow=
outdir=`pwd`/mmq$mmq

RESOURCES=~/balance/resources
GATK=$RESOURCES/GATK3.8.jar
HG19=$RESOURCES/human_g1k_v37_decoy.fasta
DBSNP=$RESOURCES/dbsnp_147_b37_common_all_20160601.vcf
BAMS=$(ls /home/ljl11/ndata1/genotyper1/science_neurons/1465/mda_*/wg.bam /home/ljl11/ndata1/genotyper1/science_neurons/1465/gatk/*.bam)

module load java

region=$(awk "NR == $SLURM_ARRAY_TASK_ID" regions.${slow}txt) \

echo $region

java -Xms9G -Xmx9G -jar $GATK \
        -mmq $mmq \
        -T HaplotypeCaller -R $HG19 \
        $(for b in $BAMS; do echo -I $b; done) \
        --dontUseSoftClippedBases \
        --dbsnp $DBSNP \
        -l INFO \
        -rf BadCigar \
        -L $region \
        -o $outdir/hc_raw.${SLURM_ARRAY_TASK_ID}.mmq${mmq}.${slow}vcf
