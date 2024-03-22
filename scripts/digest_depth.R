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
            snakemake@input['config_yaml'],
            snakemake@input['tab'],
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
if (length(args) < 4) {
    stop("usage: digest_depth.R single.cell.id config.yaml joint_depth_matrix.tab.gz out.rda [n.cores]")
}

sc.sample <- args[1]
config.yaml <- args[2]
matrix.path <- args[3]
out.rda <- args[4]

n.cores <- 1
if (length(args) == 5)
    n.cores <- as.integer(args[5])

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))

suppressMessages(library(scan2))
suppressMessages(library(future))
suppressMessages(library(progressr))
plan(multicore, workers=n.cores)

object <- make.scan(config.path=config.yaml, single.cell=sc.sample)

progressr::with_progress({
    # handler_newline causes alot of printing, but it's log-friendly
    progressr::handlers(progressr::handler_newline())
    results <- digest.depth.profile(object=object, matrix.path=matrix.path)
}, enable=TRUE)

dptab <- results$dptab
dptabs.sex <- results$dptabs.sex
clamp.dp <- results$clamp.dp
save(dptab, dptabs.sex, clamp.dp, file=out.rda)

if ('snakemake' %in% ls()) {
    sink()
}
