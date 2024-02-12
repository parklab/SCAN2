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
            snakemake@input['joint_depth_matrix'],
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
    stop("usage: spatial_sens_depth.R scan2_object.rda joint_depth_matrix.tab.gz output.txt n_cores")
}

scan2.path <- args[1]
joint.dp.path <- args[2]
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
bulk.id <- copy(object@bulk)
static.filter.params <- copy(object@static.filter.params)
genome.string <- copy(object@genome.string)
rm(object)
rm(results)
gc()
ls()

with_progress({
    handlers(handler_newline())
    depth <- compute.spatial.sensitivity.depth(single.cell.id=single.cell.id, bulk.id=bulk.id,
        static.filter.params=static.filter.params, joint.dptab.path=joint.dp.path,
        sens.tilewidth=1e3, genome.string=genome.string)
}, enable=TRUE)

cat("Writing results to", out.txt, "\n")
fwrite(depth, file=out.txt)

if ('snakemake' %in% ls()) {
    sink()
}
