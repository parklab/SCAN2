#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -c 1
#SBATCH --mem-per-cpu=4G
#SBATCH -t 120:00:00


if [ $# -ne 5 ]; then
    echo "usage: $0 working_dir single_cell_ID single_cell_bam bulk_ID bulk_bam"
    echo "Both BAMs must be indexed."
    exit 1
fi

outdir=$1
scbam=$(realpath $3)
bulkbam=$(realpath $5)

checkenv=`which check_env.sh`
if [ -z $check_env.sh ]; then
    echo "ERROR: check_env.sh script not found. The scan-snv/scripts directory"
    echo "is probably not in \$PATH."
    exit 1
fi

check_env.sh
if [ $? -ne 0 ]; then
    echo "ERROR: utils/check_env.sh failed. Please rectify all ERROR lines above."
    exit 1
else
    echo
    echo "All dependencies found. Proceeding."
    echo "-------------------------------------------------------------------------------"
fi

mkdir -p $outdir
cd $outdir

# set up links to the BAMs
check_and_link_bam()
{
    if [ $# -ne 2 ]; then
        echo "ERROR: check_and_link_bam requires 2 arguments"
        exit 1
    fi

    bam=$1
    linkname=$2

    if [ ! -f $bam ]; then
        echo "ERROR: BAM $bam not found"
        exit 1
    fi
    if [ ! -f $bam.bai ] && [ -f ${bam/.bam/.bai} ]; then
        echo "ERROR: BAM index for $bam not found"
        exit 1
    fi
    ln -s $bam $linkname
    if [ -f $bam.bai ]; then
        ln -s $bam.bai $linkname.bai
    else
        ln -s ${bam/.bam/.bai} $linkname.bai
    fi
}

pipeline_alert()
{
    echo "------------------------------------------------------------------------------"
    echo "PIPELINE: $@"
    echo "------------------------------------------------------------------------------"
}

check_and_link_bam $scbam wg.bam
check_and_link_bam $bulkbam bulk.bam

exit 1

pipeline_alert "Running GATK."
mkdir -p gatk
srun run_gatk_slurm.sh gatk wg.bam bulk.bam
if [ $? -ne 0 ]; then exit 1; fi



pipeline_alert "Running SHAPEIT2."
run_shapeit.sh hc_raw.mmq60.vcf . hunamp 22
if [ $? -ne 0 ]; then exit 1; fi

pipeline_alert "Building training sites for AB model."
get_hsnps_singlecell.sh h25 hc_raw.mmq60.vcf phased_hsnps.chr22.vcf .
if [ $? -ne 0 ]; then exit 1; fi

pipeline_alert "Fitting AB model correlation function."
ppg=$(echo 2000/$ncores | bc)
echo "DEMO: computing $ppg points per grid over $ncores cores"
mkdir -p gridfit/chr22
gridfit_chr.py \
    --bindata=training_chr22.bin \
    --local \
    --ngrids $ncores \
    --points-per-grid $ppg \
    --resume \
    --laplace=laplace_cpu_gcc \
    --outprefix=gridfit/chr22
if [ $? -ne 0 ]; then exit 1; fi

make_fits.R gridfit fits.rda
if [ $? -ne 0 ]; then exit 1; fi

pipeline_alert "Running SCAN-SNV"
mkdir -p scan-snv
scan_snv.sh hc_raw.mmq60.vcf hc_raw.mmq1.vcf . h25 hunamp scan-snv 0.1 500
if [ $? -ne 0 ]; then exit 1; fi
