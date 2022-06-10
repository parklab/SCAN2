#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
print(args)
if (length(args) != 8)
    stop('usage: call_mutations.R sc.sample config.yaml integrated_table.tab.gz abmodel_fits.rda single_cell_cigars.tab.gz bulk_cigars.tab.gz cigar_training_data.rda output.rda')

sc.sample <- args[1]
config.yaml <- args[2]
int.tab <- args[3]
abmodel.fits <- args[4]
sccigars <- args[5]
bulkcigars <- args[6]
cigardata <- args[7]
out.rda <- args[8]

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))

library(yaml)
y <- yaml::read_yaml(config.yaml)
bulk.sample <- y$bulk_sample
genome.string <- y$genome
legacy <- y$mimic_legacy

library(scan2)
library(future)
library(progressr)
plan(multicore)

with_progress({
    # handler_newline causes alot of printing, but it's log-friendly
    handlers(handler_newline())
    results <- run.pipeline(
        sc.sample=sc.sample, bulk.sample=bulk.sample,
        int.tab=int.tab, abfits=abmodel.fits,
        sccigars=sccigars, bulkcigars=bulkcigars,
        trainingcigars=cigardata,
        config.yaml=config.yaml,
        genome=genome.string,
        verbose=FALSE, report.mem=FALSE, legacy=legacy)
}, enable=TRUE)

save(results, file=out.rda, compress=FALSE)
