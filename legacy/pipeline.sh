#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -c 1
#SBATCH --mem-per-cpu=1G
#SBATCH -t 120:00:00

if [ $# -ne 5 ]; then
    echo "usage: $0 working_dir single_cell_ID single_cell_bam bulk_ID bulk_bam"
    echo "Both BAMs must be indexed."
    exit 1
fi

# Some parameters that should really be specified on the command line.
hsnp_samples=20000
laplace_chunksize=250
somatic_fdr=0.1
# setting chrs is NOT supported because most R scripts/GATK/etc assume
# all autosomes.
# when it works: chrs must be a SPACE SEPARATED list of chromosomes
# preferably in increasing order.
# chrX will one day be supported as chr=23 because we use SLURM array
# indeces, which must be integers, to parallelize over chroms
chrs=`seq 1 22`
arraychrs=`echo $chrs | tr ' ' ','`

# Exit if any individual command fails. Report the last command on exit
set -e
trap 'echo "Failed command (exit=$?): $BASH_COMMAND"' EXIT

outdir=$1
scid=$2
scbam=$(realpath $3)
bulkid=$4
bulkbam=$(realpath $5)

pipeline_alert()
{
    echo "-------------------------------------------------------------------------------"
    echo "PIPELINE: $@"
    echo "-------------------------------------------------------------------------------"
}


pipeline_alert "Checking environment variables."
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

check_and_link_bam $scbam wg.bam
check_and_link_bam $bulkbam bulk.bam


pipeline_alert "Running GATK."
mkdir -p gatk
#sbatch -W -o gatk/log.%a.txt run_gatk_slurm.sh gatk wg.bam bulk.bam

# Don't use globbing; need to retain 1-200 order
#sbatch -W -o gatk/cat60.log cat_vcfs.sh \
#    gatk/hc_raw.mmq60.vcf \
#    $(for i in `seq 1 205`; do echo "gatk/hc_raw.$i.mmq60.vcf"; done)
#sbatch -W -o gatk/cat1.log cat_vcfs.sh \
#    gatk/hc_raw.mmq1.vcf \
#    $(for i in `seq 1 205`; do echo "gatk/hc_raw.$i.mmq1.vcf"; done)


pipeline_alert "Running SHAPEIT2."
mkdir -p shapeit
sbatch --array=$arraychrs -W -o shapeit/log.%a.txt \
    run_shapeit2.sh gatk/hc_raw.mmq60.vcf shapeit $bulkid
sbatch -W -o shapeit/cathsnps.log \
    cat_vcfs.sh phased_hsnps.vcf \
        $(for i in $chrs; do echo shapeit/phased_hsnps.chr$i.vcf; done)


pipeline_alert "Building training sites for AB model."
mkdir -p ab_model
# TODO: get_hsnps_singlecell.sh calls torda.R, which creates training data
# for chrs 1-22 despite user settings. Won't create a file for chrX.
sbatch -W -o get_hsnps.log \
    get_hsnps_singlecell.sh $scid \
        gatk/hc_raw.mmq60.vcf phased_hsnps.vcf ab_model


pipeline_alert "Fitting AB model correlation functions for each chromosome."
sbatch -W -o ab_model/log.txt \
    gridfit_slurm.sh ab_model $laplace_chunksize park $chrs

# Check all the output files
for chrom in $chrs; do
    if [ ! -f "ab_model/gridfit/chr$chrom/fit.rda" ]; then
        echo "Failed (at least) chromosome $chrom"
        exit 1
    fi
done

# TODO: make_fits.R uses chr1-22 static; cannot handle chrX
make_fits.R ab_model/gridfit ab_model/fits.rda



pipeline_alert "Running SCAN-SNV"
mkdir -p scan-snv
sbatch -W -o scan-snv/log.txt \
    scan_snv.sh \
    gatk/hc_raw.mmq60.vcf \
    gatk/hc_raw.mmq1.vcf \
    . \
    $scid \
    $bulkid \
    scan-snv \
    $somatic_fdr \
    $hsnp_samples \
    ab_model
