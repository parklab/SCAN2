#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
print(args)
if (length(args) < 10)
    stop('usage: run_pipeline.R sc.sample bulk.sample mmq60.tab.gz mmq1.tab.gz hsnps.tab.gz abmodel_fits.rda single_cell_cigars.tab.gz bulk_cigars.tab.gz genome_string output.rda [tmpsave.rda]')

sc.sample <- args[1]
bulk.sample <- args[2]
mmq60 <- args[3]
mmq1 <- args[4]
hsnps <- args[5]
abmodel.fits <- args[6]
sccigars <- args[7]
bulkcigars <- args[8]
genome.string <- args[9]
out.rda <- args[10]

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))


library(scan2)
library(future)
library(progressr)
plan(multicore)

if (length(args) == 11) {
    tmpsave.rda <- args[11]
with_progress(
    results <- run.pipeline(
        sc.sample=sc.sample, bulk.sample=bulk.sample,
        mmq60=mmq60, mmq1=mmq1,
        hsnps=hsnps, abfits=abmodel.fits,
        sccigars=sccigars, bulkcigars=bulkcigars,
        tmpsave.rda=tmpsave.rda, genome=genome.string, verbose=TRUE)
, enable=TRUE)
} else {
with_progress(
    results <- run.pipeline(
        sc.sample=sc.sample, bulk.sample=bulk.sample,
        mmq60=mmq60, mmq1=mmq1,
        hsnps=hsnps, abfits=abmodel.fits,
        sccigars=sccigars, bulkcigars=bulkcigars,
        genome=genome.string, verbose=TRUE)
, enable=TRUE)
}

save(results, file=out.rda, compress=FALSE)
