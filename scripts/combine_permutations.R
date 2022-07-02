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
            snakemake@input['permrdas'],
            snakemake@params['genome_string'],
            snakemake@output['perm'],
            snakemake@output['seedinfo'],
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
    stop("usage: combine_permutations.R genome_string perms_output.rda seedinfo_output.rda n.cores perms1.rda [ perms2.rda ... permsN.rda ]")
}

genome.string <- args[1]
out.rda <- args[2]
out.seedinfo <- args[3]
n.cores <- as.integer(args[4])
in.rdas <- args[-(1:4)]

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))

if (file.exists(out.seedinfo))
    stop(paste('output file', out.seedinfo, 'already exists, please delete it first'))


suppressMessages(library(scan2))
suppressMessages(library(future))
suppressMessages(library(progressr))
plan(multicore, workers=n.cores)

progressr::with_progress({
    # handler_newline causes alot of printing, but it's log-friendly
    progressr::handlers(progressr::handler_newline())
    permdata <- combine.permutations(perm.files=in.rdas,
        genome=genome.string, report.mem=TRUE)
}, enable=TRUE)

zperml <- permdata$zperml
save(zperml, file=out.rda, compress=FALSE)
seed.info <- permdata$seed.info
save(seed.info, file=out.seedinfo, compress=FALSE)

if ('snakemake' %in% ls()) {
    sink()
}
