#!/bin/bash

which laplace_cpu_gcc >& /dev/null
if [ ! $? ]; then
    echo "DEMO: could not find laplace_cpu_gcc in PATH, please compile the gcc"
    echo "      version of the Laplace approximator."
    exit 1
fi

if [ -z ${PATH+x} ]; then
    echo "DEMO: please set the \$PATH environment variable"
    exit 1
fi
if [ -z ${LD_LIBRARY_PATH+x} ]; then
    echo "DEMO: please set the \$LD_LIBRARY_PATH environment variable"
    exit 1
fi
if [ -z ${GATK_PATH+x} ]; then
    echo "DEMO: please set the \$GATK_PATH environment variable"
    exit 1
fi
if [ -z ${SHAPEIT_ROOT+x} ]; then
    echo "DEMO: please set the \$SHAPEIT_ROOT environment variable"
    exit 1
fi
if [ -z ${REFPANEL_ROOT+x} ]; then
    echo "DEMO: please set the \$REFPANEL_ROOT environment variable"
    exit 1
fi

ncores=8

mkdir demo
cd demo

echo "DEMO: downloading read data.. (~2 GB)"
wget http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam
wget http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam.bai
wget http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam
wget http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam.bai
ln -s il-12.chr22.bam wg.bam
ln -s il-12.chr22.bam.bai wg.bam.bai

echo "DEMO: running GATK.. (~1.5 hours with 8 threads)"
run_gatk_demo.sh $ncores 60 . hunamp.chr22.bam il-12.chr22.bam
run_gatk_demo.sh $ncores 1 . hunamp.chr22.bam il-12.chr22.bam

echo "DEMO: phasing hSNPs with SHAPEIT2"
run_shapeit.sh hc_raw.mmq60.vcf . hunamp 22

echo "DEMO: building training sites for AB model.."
get_hsnps_singlecell.sh h25 hc_raw.mmq60.vcf phased_hsnps.chr22.vcf .

echo 'DEMO: fitting AB model correlation function.. (~35 minutes with 8 threads)'
mkdir -p gridfit/chr22
gridfit_chr.py \
    --bindata=training_chr22.bin \
    --local \
    --ngrids $ncores \
    --points-per-grid 250 \
    --resume \
    --laplace=laplace_cpu_gcc \
    --outprefix=gridfit/chr22

make_fits.R gridfit fits.rda

echo "DEMO: running SCAN-SNV"
mkdir -p scan-snv
scan_snv.sh hc_raw.mmq60.vcf hc_raw.mmq1.vcf . h25 hunamp scan-snv 0.1
