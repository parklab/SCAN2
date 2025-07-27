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
            snakemake@params['genome'],
            snakemake@output['rda'],
            snakemake@input[['tabs']]
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("usage: summarize_depth_gather.R genome out.rda in_depth_table1.txt [ in_depth_table2.txt ... in_depth_tableN.txt ]")
}

genome <- args[1]
out.rda <- args[2]
in.txts <- args[-(1:2)]

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))

suppressMessages(library(scan2))

sex.chroms <- genome.string.to.chroms(genome, group='sex')

dptabs <- lapply(in.txts, function(fname) {
    print(fname)
    header <- strsplit(strsplit(readLines(fname, n=1), '\t')[[1]], '=')
    vals <- as.list(setNames(sapply(header, `[`, 2), sapply(header, `[`, 1)))
    m <- as.matrix(read.csv(fname, skip=1, header=FALSE, sep='\t'))
    dimnames(m) <- setNames(list(0:vals$sc_max_depth, 0:vals$bulk_max_depth),
        c(vals$sc_sample, vals$bulk_sample))
    vals$dptab <- as.table(m)
    vals
})

names(dptabs) <- sapply(dptabs, function(dpt) dpt$chrom)

dptab <- Reduce(`+`, lapply(dptabs[!(names(dptabs) %in% sex.chroms)], function(dpt) dpt$dptab))
dptabs.sex <- lapply(dptabs[names(dptabs) %in% sex.chroms], function(dpt) dpt$dptab)
# Currently there's no support for different max_depth for single cell and bulk.
# Also, it is also expected to be uniform across chromosomes. So just take the first one.
clamp.dp <- as.integer(dptabs[[1]]$sc_max_depth)  # the header k=v splitting doesn't detect types
save(dptab, dptabs.sex, clamp.dp, file=out.rda)

if ('snakemake' %in% ls()) {
    sink()
}
