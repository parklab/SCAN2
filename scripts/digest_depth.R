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
            snakemake@input['tab'],
            snakemake@params['sc_sample'],
            snakemake@params['bulk_sample'],
            snakemake@params['genome'],
            snakemake@output['rda'],
            snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
print(args)
if (length(args) < 5) {
    stop("usage: digest_depth.R joint_depth_matrix.tab.gz sc.sample bulk.sample genome out.rda [n.cores]")
}

path <- args[1]
sc.sample <- args[2]
bulk.sample <- args[3]
genome <- args[4]
out.rda <- args[5]

n.cores <- 1
if (length(args) == 6)
    n.cores <- as.integer(args[6])

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))

suppressMessages(library(scan2))
suppressMessages(library(future))
suppressMessages(library(progressr))
plan(multicore, workers=n.cores)

progressr::with_progress({
    # handler_newline causes alot of printing, but it's log-friendly
    progressr::handlers(progressr::handler_newline())
    results <- digest.depth.profile(path=path, sc.sample=sc.sample, bulk.sample=bulk.sample,
        genome=genome)
}, enable=TRUE)

dptab <- results@dptab
clamp.dp <- results@clamp.dp
save(dptab, clamp.dp, file=out.rda)

if ('snakemake' %in% ls()) {
    sink()
}
