#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
print(args)
if (length(args) < 16) {
    cat("To compute AB estimates, models and excess cigar scores as part of this job, supply the strings 'NULL' as the file paths.\n")
    cat("For each of the final optional arguments (n.cores, scan2.version...), all preceding optional args must be specified.\n")
    stop('usage: call_mutations.R sc.sample config.yaml integrated_table.tab.gz abmodel_fits.rda single_cell_cigars.tab.gz bulk_cigars.tab.gz cigar_training_data.rda dptab.rda sc.binned.counts bulk.binned.counts gc.bins precomputed.ab.and.models.tab.gz precomputed.excess.cigar.scores.tab.gz abmodel_covariates.tab.gz depth_covariates.tab.gz out.rda [n.cores] [scan2.version] [scan2.buildnum] [scan2.githash]')
}

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
precomputed.ab.and.models <- ifelse(args[12] == "NULL", NULL, args[12])
precomputed.excess.cigar.scores <- ifelse(args[13] == "NULL", NULL, args[13])
abmodel.covs <- args[14]
depth.covs <- args[15]
out.rda <- args[16]

n.cores <- 1
if (length(args) >= 17)
    n.cores <- as.integer(args[17])

scan2.version <- NA
if (length(args) >= 18)
    scan2.version <- args[18]

scan2.buildnum <- NA
if (length(args) >= 19)
    scan2.buildnum <- args[19]

scan2.githash <- NA
if (length(args) == 20)
    scan2.githash <- args[20]

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))

library(scan2)
library(future)
library(progressr)
plan(multicore, workers=n.cores)

results <- make.scan(config.path=config.yaml,
        single.cell=sc.sample,
        # Record version info
        pipeline.version=scan2.version,
        pipeline.buildnum=scan2.buildnum,
        pipeline.githash=scan2.githash)

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
        ab.ests.and.models.path=precomputed.ab.and.models,
        excess.cigar.scores.path=precomputed.excess.cigar.scores,
        abmodel.covs=abmodel.covs,
        depth.covs=depth.covs,
        # want coarse tiles for this job.
        grs.for.parallelization=analysis.set.tiling.for.parallelization(results, total.tiles=24),
        verbose=FALSE, report.mem=TRUE)
}, enable=TRUE)

save(results, file=out.rda, compress=FALSE)
