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
            snakemake@params['min_sc_dp'],
            snakemake@params['min_bulk_dp'],
            snakemake@output['bed'],
            snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
print(args)
if (length(args) < 7) {
    stop("usage: make_callable_regions_bed.R joint_depth_matrix.tab.gz sc.sample bulk.sample genome min_sc_dp min_bulk_dp out.bed [n.cores]")
}

path <- args[1]
sc.sample <- args[2]
bulk.sample <- args[3]
genome <- args[4]
min.sc.dp <- as.integer(args[5])
min.bulk.dp <- as.integer(args[6])
out.bed <- args[7]

n.cores <- 1
if (length(args) == 8)
    n.cores <- as.integer(args[8])

if (file.exists(out.bed))
    stop(paste('output file', out.bed, 'already exists, please delete it first'))

suppressMessages(library(scan2))
suppressMessages(library(future))
suppressMessages(library(progressr))
plan(multicore, workers=n.cores)

progressr::with_progress({
    # handler_newline causes alot of printing, but it's log-friendly
    progressr::handlers(progressr::handler_newline())
    results <- make.callable.regions(path=path, sc.sample=sc.sample, bulk.sample=bulk.sample,
        genome=genome, min.sc.dp=min.sc.dp, min.bulk.dp=min.bulk.dp)
}, enable=TRUE)


write.table(head(as.data.frame(results)[,1:3]),
    file=out.bed, row.names=F, col.names=F, quote=F, sep='\t')

if ('snakemake' %in% ls()) {
    sink()
}
