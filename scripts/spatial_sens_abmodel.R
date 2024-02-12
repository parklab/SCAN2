#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) {
        ret <- unlist(c(
            snakemake@input['scan2_object'],
            snakemake@input['integrated_table'],
            snakemake@output['txt'],
            snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
    stop("usage: spatial_sens_abmodel.R scan2_object.rda integrated_table.tab.gz output.txt n_cores")
}

scan2.path <- args[1]
integrated.table.path <- args[2]
out.txt <- args[3]
n.cores <- as.integer(args[4])

for (f in out.txt)
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))

suppressMessages(library(scan2))
suppressMessages(library(future))
suppressMessages(library(progressr))
plan(multicore, workers=n.cores)


# MUST delete this object before running the pipeline to prevent copying the
# memory to child processes.
object <- get(load(scan2.path)) # loads "results"
single.cell.id <- copy(object@single.cell)
ab.fits <- copy(object@ab.fits)
genome.string <- copy(object@genome.string)
rm(object)
rm(results)
gc()
ls()

with_progress({
    handlers(handler_newline())
    abmodel <- compute.spatial.sensitivity.abmodel(single.cell.id=single.cell.id,
        ab.fits=ab.fits, integrated.table.path=integrated.table.path,
        sens.tilewidth=1e3, genome.string=genome.string)
}, enable=TRUE)

cat("Writing results to", out.txt, "\n")
fwrite(abmodel, file=out.txt)

if ('snakemake' %in% ls()) {
    sink()
}
