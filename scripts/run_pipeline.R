#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
print(args)
if (length(args) < 9)
    stop('usage: run_pipeline.R sc.sample bulk.sample integrated_table.tab.gz abmodel_fits.rda single_cell_cigars.tab.gz bulk_cigars.tab.gz cigar_training_data.rda genome_string output.rda [tmpsave.rda]')

sc.sample <- args[1]
bulk.sample <- args[2]
int.tab <- args[3]
abmodel.fits <- args[4]
sccigars <- args[5]
bulkcigars <- args[6]
cigardata <- args[7]
genome.string <- args[8]
out.rda <- args[9]

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))


library(scan2)
library(future)
library(progressr)
plan(multicore)

if (length(args) == 10) {
    tmpsave.rda <- args[10]
    with_progress({
        handlers(handler_newline())
        results <- run.pipeline(
            sc.sample=sc.sample, bulk.sample=bulk.sample,
            int.tab=int.tab, abfits=abmodel.fits,
            sccigars=sccigars, bulkcigars=bulkcigars,
            trainingcigars=cigardata, fdr.prior.data=fdrpriordata,
            tmpsave.rda=tmpsave.rda, genome=genome.string, verbose=FALSE, report.mem=FALSE, legacy=FALSE)
    }, enable=TRUE)
} else {
    with_progress({
        handlers(handler_newline())
        results <- run.pipeline(
            sc.sample=sc.sample, bulk.sample=bulk.sample,
            int.tab=int.tab, abfits=abmodel.fits,
            sccigars=sccigars, bulkcigars=bulkcigars,
            trainingcigars=cigardata, fdr.prior.data=fdrpriordata,
            genome=genome.string, verbose=FALSE, report.mem=FALSE, legacy=FALSE)
    }, enable=TRUE)
}

save(results, file=out.rda, compress=FALSE)
