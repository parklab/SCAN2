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
            fits=snakemake@output['fits'],
            fit.details=snakemake@output['fit_details'],
            snakemake@input   # one or more
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("usage: abmodel_gather_tiled_script.R out_fits.rda out_details.rda in_fit1.rda [ in_fit2.rda .. in_fitN.rda ]")
}

out.fits <- args[1]
out.details <- args[2]
in.rdas <- args[-(1:2)]

if (file.exists(out.fits))
    stop(paste('output file', out.fits, 'already exists, please delete it first'))

if (file.exists(out.details))
    stop(paste('output file', out.details, 'already exists, please delete it first'))

# Each load() produces a list with one element named after the chromosome fit
details <- do.call(c, lapply(in.rdas, function(f) get(load(f))))
fits <- do.call(rbind, lapply(details, function(d) d$fit))

save(fits, file=out.fits)
save(details, file=out.details)

if ('snakemake' %in% ls()) {
    sink()
}
