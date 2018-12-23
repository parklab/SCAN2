#!/bin/bash
#SBATCH -c 2
#SBATCH --mem-per-cpu=2G
#SBATCH -p medium
#SBATCH -t 120:00:00


# just an example use of runchr.py

if [ $# -ne 1 -a $# -ne 2 ] ; then
    echo "usage: $0 dir [queue]"
    echo "dir must contain binary training data in the proper format, named chrXXX.bin"
    exit 1
fi

dir=$1
queue=park
if [ $# -eq 2 ]; then
    queue=$2
fi


module load gcc
module load openblas
module load slurm-drmaa
module load R/3.3.3

export DRMAA_LIBRARY_PATH=/n/app/slurm-drmaa/1.0.7/lib/libdrmaa.so

for chrom in `seq 1 22`; do
    f="$dir/chr$chrom.bin"
    if [ ! -f $f ]; then
        echo "$0: expected binary file $f not present"
        exit 1
    fi
done

for chrom in `seq 1 22`; do
    # resume mode is always safe to use, even for new analyses
    cmd="gridfit_chr.py
        --bindata=$dir/chr$chrom.bin \
        --outprefix=$dir/gridfit/chr$chrom/ \
        --laplace=laplace_cpu_gcc \
        --queue=$queue \
        --resume"
    mkdir -p $dir/gridfit/chr$chrom
    echo $cmd
    stdbuf -o0 -e0 $cmd &>> $dir/gridfit/chr$chrom/log.txt &
done

# waits for all child processes
wait

for chrom in `seq 1 22`; do
    if [ ! -f "$dir/gridfit/chr$chrom/fit.rda" ]; then
        echo "Failed (at least) chromosome $chrom"
        exit 1
    fi
done

# If we got here, then should be good
make_fits.R $dir/gridfit $dir/fits.rda
