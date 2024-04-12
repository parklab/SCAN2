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
            snakemake@input['rda'],
            snakemake@output['rda'],
            snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("usage: summarize_object.R scan2_object.rda out.rda [ n.cores ]")
}

in.rda <- args[1]
out.rda <- args[2]

n.cores <- 1
if (length(args) > 2)
    n.cores <- as.integer(args[3])

for (f in c(out.rda)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(future))
suppressMessages(library(progressr))

if (n.cores > 1) {
    plan(multicore, workers=n.cores)
} else {
    plan(sequential)
}
cat("Summarizing SCAN2 object with", nbrOfWorkers(), "cores.\n")

progressr::handlers(global=TRUE)

load(in.rda, verb=TRUE)
results@static.filter.params$mode <- 'new'  # XXX: deleteme: old bug workaround
smry <- make.summary.scan2(results, quiet=FALSE, preserve.object=FALSE)
save(smry, file=out.rda, compress=FALSE)

if ('snakemake' %in% ls()) {
    sink()
}
