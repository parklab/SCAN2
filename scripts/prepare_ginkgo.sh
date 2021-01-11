#!/bin/bash

if [ $# -ne 2 ]; then
    echo "usage: $0 scan2_dir output_dir"
    exit 1
fi

ssdir=$1
outdir=$2

ssconfig=$ssdir/scan.yaml


scriptdir=$(cat $ssconfig | shyaml get-value scripts)
bamtobedscript=$scriptdir/run_bamtobed.sh
ginkgoscript=$scriptdir/run_ginkgo.sh

if [ ! -f $ssconfig ]; then
    echo "ERROR: SCAN2 directory $ssdir does not contain a 'scan.yaml'"
    exit 1
fi

echo "WARNING: be sure to 'conda activate ginkgo' before running this script!"

mkdir -p $outdir 

bams=$(cat $ssconfig | shyaml values bam_map)

tmpdir=`pwd`
cd $outdir
for bam in $bams; do
    sbatch $bamtobedscript $bam $(basename $bam .bam).bed
done
cd $tmpdir
