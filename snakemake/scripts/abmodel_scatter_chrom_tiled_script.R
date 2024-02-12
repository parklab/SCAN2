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
            snakemake@input['inttab'],
            snakemake@input['config'],
            snakemake@params['sc_sample'],
            snakemake@params['chrom'],
            snakemake@output['rda'],
            snakemake@params['n_tiles'],
            n_cores=snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
    cat('if n.cores is unspecified, default is 1\n')
    stop("usage: abmodel_scatter_chrom_tiled_script.R integrated_table.tab.bgz config.yaml sc_sample_id chromosome out.rda n.tiles [n.cores]")
}

int.tab.bgz <- args[1]
config.yaml <- args[2]
sc.sample <- args[3]
chrom <- args[4]
out.rda <- args[5]
n.tiles <- as.integer(args[6])
library(future)
if (length(args) == 7) {
    plan(multicore)
    n.cores <- as.integer(args[7])
} else {
    plan(sequential)
    n.cores <- 1
}

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))


suppressMessages(library(scan2))
suppressMessages(library(future))

cat(future::availableCores(), "cores detected by library(future),", n.cores, "requested by user", "\n")
plan(multicore, workers=n.cores)

s <- make.scan(config.path=config.yaml, single.cell=sc.sample)

f <- compute.ab.fits(s, path=int.tab.bgz, chroms=chrom, n.tiles=n.tiles)

save(f, file=out.rda)

if ('snakemake' %in% ls()) {
    sink()
}
