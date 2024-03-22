#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
print(args)
if (length(args) < 12)
    stop('usage: call_mutations.R sc.sample config.yaml integrated_table.tab.gz abmodel_fits.rda single_cell_cigars.tab.gz bulk_cigars.tab.gz cigar_training_data.rda dptab.rda sc.binned.counts bulk.binned.counts gc.bins out.rda [n.cores]')

sc.sample <- args[1]
config.yaml <- args[2]
int.tab <- args[3]
abmodel.fits <- args[4]
sccigars <- args[5]
bulkcigars <- args[6]
cigardata <- args[7]
dptab <- args[8]
sc.bins <- args[9]
bulk.bins <- args[10]
gc.bins <- args[11]
out.rda <- args[12]

n.cores <- 1
if (length(args) == 13)
    n.cores <- as.integer(args[13])

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))


library(scan2)
library(future)
library(progressr)
plan(multicore, workers=n.cores)

results <- make.scan(config.path=config.yaml, single.cell=sc.sample)

with_progress({
    # handler_newline causes alot of printing, but it's log-friendly
    handlers(handler_newline())
    results <- run.pipeline(results,
        int.tab=int.tab, abfits=abmodel.fits,
        sccigars=sccigars, bulkcigars=bulkcigars,
        trainingcigars=cigardata, dptab=dptab,
        sc.binned.counts=sc.bins,
        bulk.binned.counts=bulk.bins,
        gc.content.bins=gc.bins,
        verbose=FALSE, report.mem=TRUE)
}, enable=TRUE)

save(results, file=out.rda, compress=FALSE)
