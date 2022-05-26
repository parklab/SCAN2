#!/usr/bin/env Rscript


# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) {
stop('this script has not been adapted for snakemake yet')
        ret <- unlist(c(
            snakemake@input['hsnps'],
            snakemake@params['chrom'],
            snakemake@output['rda']
            snakemake@params['n_tiles'],
            snakemake@params['n_cores']
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    cat('if n.cores is unspecified, default is 1\n')
    stop("usage: abmodel_scatter_chrom_tiled_script.R hsnps.tab.bgz chromosome out.rda n.tiles [n.cores]")
}

hsnps.tab.bgz <- args[1]
chrom <- args[2]
out.rda <- args[3]
n.tiles <- as.integer(args[4])
if (length(args) == 5) {
    library(future)
    plan(multicore)
    n.cores <- as.integer(args[5])
} else {
    plan(sequential)
    n.cores <- 1
}

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))


suppressMessages(library(scan2))

s <- make.scan(bulk='dummy.bulk', genome='hs37d5', single.cell='dummy.sc')

f <- compute.ab.fits(s, path=hsnps.tab.bgz, chroms=chrom, n.tiles=n.tiles)

save(f, file=out.rda)

if ('snakemake' %in% ls()) {
    sink()
}
