#!/bin/bash

if [ $# -ne 5 ]; then
    echo "usage: $0 gatk.tab.gz hsnps.tab.gz scan2_config.yaml out.tab.gz out_details.rda"
    exit 1
fi

gatk=$1
hsnps=$2
scan2config=$3
outtab=$4
outrda=$5

if [ -f $outtab ]; then
    echo "output file $outtab already exists, please delete it first"
    exit 1
fi

if [ -f $outrda ]; then
    echo "output file $outrda already exists, please delete it first"
    exit 1
fi


bulk_id=$(cat $scan2config | shyaml get-value bulk_sample)
sc_id=$(cat $scan2config  | shyaml keys sc_bams | head -1) # Get first sample. doesn't matter which as long as it's a valid sample

# This reads GATK, the hSNPs tab, and the SCAN2 config file and uses
# SCAN2 methods to downsample hSNPs. The resulting table is written to
# stdout for bgzip and tabix.
# NOTE: the genome type in make.scan doesn't matter at all
Rscript -e "suppressMessages(library(scan2)); s <- make.scan('$sc_id', '$bulk_id', 'hs37d5'); s <- add.static.filter.params(s, '$scan2config'); s <- read.gatk(s, '$gatk', quiet=TRUE, add.mutsig=FALSE); s <- add.training.data(s, '$hsnps', quiet=TRUE); s <- resample.training.data(s); resample.data <- s@resampled.training.data; save(resample.data, file='$outrda'); out <- s@gatk[training.site==TRUE,.(chr, pos, refnt, altnt, training.hap1, training.hap2, training.phgt, resampled.training.site)]; colnames(out) <- c('#chr', 'pos', 'refnt', 'altnt', 'hap1', 'hap2', 'phgt', 'resampled'); fwrite(out, file='/dev/stdout', sep='\t')" | bgzip -c /dev/stdin > $outtab

tabix -p vcf -S 1 $outtab
