#!/bin/bash

if [ $# -ne 2 ]; then
    echo "usage: $0 in.rda out.rda"
    exit 1
fi

inrda=$1
outrda=$2

if [ -f $outrda ]; then
    echo "output file $outrda already exists, please delete it first"
    exit 1
fi

Rscript -e 'library(scan2); load("'$inrda'"); results@gatk <- results@gatk[(pass == TRUE | rescue == TRUE) & somatic.candidate==TRUE]; save(results, file="'$outrda'")'
