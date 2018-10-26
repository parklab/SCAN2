#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: $0 n_threads"
    echo
    echo "Please specify the number of threads to use. E.g.,"
    echo "$ demo.sh 8"
    echo "to run with 8 threads. Runtime is approximately 5 hours with 8 threads"
    exit 1
fi

ncores=$1

if [ ! -f "utils/check_env.sh" ]; then
    echo "ERROR: please run the demo script from the root of the scan-snv repository"
    exit 1
fi

utils/check_env.sh
if [ $? -ne 0 ]; then
    echo "ERROR: utils/check_env.sh failed. Please rectify all ERROR lines above."
    exit 1
else
    echo
    echo "All dependencies found. Proceeding."
    echo "-------------------------------------------------------------------------------"
fi

mkdir -p demo
cd demo

echo "DEMO: downloading read data.. (~2 GB)"
if [ ! -f hunamp.chr22.bam ]; then
    wget http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam
fi
if [ ! -f hunamp.chr22.bam.bai ]; then
    wget http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam.bai
fi
if [ ! -f il-12.chr22.bam ]; then
    wget http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam
    ln -s il-12.chr22.bam wg.bam
fi
if [ ! -f il-12.chr22.bam.bai ]; then
    wget http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam.bai
    ln -s il-12.chr22.bam.bai wg.bam.bai
fi

echo "DEMO: running GATK.. (~1.5 hours with 8 threads)"
run_gatk_demo.sh $ncores 60 . hunamp.chr22.bam il-12.chr22.bam
if [ $? -ne 0 ]; then exit 1; fi
run_gatk_demo.sh $ncores 1 . hunamp.chr22.bam il-12.chr22.bam
if [ $? -ne 0 ]; then exit 1; fi

echo "DEMO: phasing hSNPs with SHAPEIT2"
run_shapeit.sh hc_raw.mmq60.vcf . hunamp 22
if [ $? -ne 0 ]; then exit 1; fi

echo "DEMO: building training sites for AB model.."
get_hsnps_singlecell.sh h25 hc_raw.mmq60.vcf phased_hsnps.chr22.vcf .
if [ $? -ne 0 ]; then exit 1; fi

echo 'DEMO: fitting AB model correlation function using 2000 samples.. (~35 minutes with 8 threads)'
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

echo "DEMO: running SCAN-SNV"
mkdir -p scan-snv
scan_snv.sh hc_raw.mmq60.vcf hc_raw.mmq1.vcf . h25 hunamp scan-snv 0.1
if [ $? -ne 0 ]; then exit 1; fi
