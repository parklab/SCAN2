#!/bin/bash
#SBATCH -c 2
#SBATCH --mem-per-cpu=4G
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 120:00:00

# This script runs 22 instances of gridfit_chr.py. Each of these scripts
# control and monitor SLURM jobs through a full grid estimate for one
# chromosome. All heavy computing is submitted by these controllers to
# SLURM; however, because 22 controller processes are being run in this
# single SLURM job, we ask for 2 cores to (hopefully) prevent being
# cancelled due to too much CPU usage.

if [ $# -lt 4 ] ; then
    echo "usage: $0 dir chunksize queue chromosome1 [ chrom2 ... chromN ]"
    echo "dir must contain binary training data in the proper format,"
    echo "named training_chrXXX.bin"
    exit 1
fi

dir=$1
chunksize=$2
queue=$3
shift 3
chrs=$@

for chrom in $chrs; do
    f="$dir/training_chr$chrom.bin"
    if [ ! -f $f ]; then
        echo "$0: expected training file $f not present"
        exit 1
    fi
done

for chrom in $chrs; do
    # resume mode is always safe to use, even for new analyses
    cmd="gridfit_chr.py
        --bindata=$dir/training_chr$chrom.bin \
        --laplace=laplace_cpu_gcc \
        --chunksize=$chunksize \
        --outprefix=$dir/gridfit/chr$chrom/ \
        --queue=$queue \
        --resume"
    mkdir -p $dir/gridfit/chr$chrom
    echo $cmd
    stdbuf -o0 -e0 $cmd &>> $dir/gridfit/chr$chrom/log.txt &
done

# waits for all child processes
wait
