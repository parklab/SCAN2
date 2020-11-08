#!/bin/bash

#SBATCH -t 32:00:00
#SBATCH --mem=8G
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -c 1

if [ $# -ne 1 ]; then
    echo "usage: $0 binsize"
    echo "only certain binsizes are allowed"
    exit 1
fi

binsize=$1

echo "WARNING!! MUST RUN 'conda activate ginkgo' FIRST!!"

~/ndata1/ginkgo/ginkgo/cli/ginkgo.sh \
    --input `pwd` \
    --genome hg19 \
    --binning \
    variable_${binsize}_150_bwa \
    --maskpsrs
