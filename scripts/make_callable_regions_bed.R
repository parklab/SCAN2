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
            snakemake@params['sc_sample'],
            snakemake@input['config'],
            snakemake@input['tab'],
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
if (length(args) < 4) {
    stop("usage: make_callable_regions_bed.R sc.sample config.yaml joint_depth_matrix.tab.gz out.bed [n.cores]")
}

sc.sample <- args[1]
muttype <- args[2]
config.path <- args[3]
matrix.path <- args[4]
out.bed <- args[5]

n.cores <- 1
if (length(args) == 6)
    n.cores <- as.integer(args[6])

if (file.exists(out.bed))
    stop(paste('output file', out.bed, 'already exists, please delete it first'))

suppressMessages(library(scan2))
suppressMessages(library(future))
suppressMessages(library(progressr))
plan(multicore, workers=n.cores)

object <- make.scan(single.cell=sc.sample, config.path=config.path)

progressr::with_progress({
    # handler_newline causes alot of printing, but it's log-friendly
    progressr::handlers(progressr::handler_newline())
    results <- make.callable.regions(object=object, matrix.path=matrix.path, muttype=muttype)
}, enable=TRUE)


write.table(as.data.frame(results)[,1:3],
    file=out.bed, row.names=F, col.names=F, quote=F, sep='\t')

if ('snakemake' %in% ls()) {
    sink()
}
