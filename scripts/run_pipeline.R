#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
print(args)
if (length(args) < 12)
    stop('usage: run_pipeline.R sc.sample bulk.sample mmq60.tab.gz mmq1.tab.gz hsnps.tab.gz abmodel_fits.rda single_cell_cigars.tab.gz bulk_cigars.tab.gz cigar_training_data.rda fdr_prior_data.rda genome_string output.rda [tmpsave.rda]')

sc.sample <- args[1]
bulk.sample <- args[2]
mmq60 <- args[3]
mmq1 <- args[4]
hsnps <- args[5]
abmodel.fits <- args[6]
sccigars <- args[7]
bulkcigars <- args[8]
cigardata <- args[9]
fdrpriordata <- args[10]
genome.string <- args[11]
out.rda <- args[12]

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))


library(scan2)
library(future)
library(progressr)
plan(multicore)

if (length(args) == 13) {
    tmpsave.rda <- args[13]
    with_progress(
        results <- run.pipeline(
            sc.sample=sc.sample, bulk.sample=bulk.sample,
            mmq60=mmq60, mmq1=mmq1,
            hsnps=hsnps, abfits=abmodel.fits,
            sccigars=sccigars, bulkcigars=bulkcigars,
            trainingcigars=cigardata, fdr.prior.data=fdrpriordata,
            tmpsave.rda=tmpsave.rda, genome=genome.string, verbose=TRUE)
    , enable=TRUE)
} else {
    with_progress(
        results <- run.pipeline(
            sc.sample=sc.sample, bulk.sample=bulk.sample,
            mmq60=mmq60, mmq1=mmq1,
            hsnps=hsnps, abfits=abmodel.fits,
            sccigars=sccigars, bulkcigars=bulkcigars,
            trainingcigars=cigardata, fdr.prior.data=fdrpriordata,
            genome=genome.string, verbose=TRUE)
    , enable=TRUE)
}

save(results, file=out.rda, compress=FALSE)
