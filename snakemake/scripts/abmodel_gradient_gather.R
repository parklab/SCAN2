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
    stop("usage: abmodel_gradient_gather.R out_fits.rda out_details.rda in_fit1.rda [ in_fit2.rda .. in_fitN.rda ]")
}

out.fits <- args[1]
out.details <- args[2]
in.rdas <- args[-(1:2)]

if (file.exists(out.fits))
    stop(paste('output file', out.fits, 'already exists, please delete it first'))

if (file.exists(out.details))
    stop(paste('output file', out.details, 'already exists, please delete it first'))

library(data.table)

# Each load() produces a list with one element named fits. fits has one element
# for each random initiation of starting parameters. since the MLE is the fit
# we want, the random starting params that produce the maximal logp value is the
# best fit.
all.data <- lapply(in.rdas, function(f) {
    load(f, verb=TRUE)  # loads 'fits', 'chrom'
    fittab <- rbindlist(lapply(fits, function(f) cbind(f[[1]], convergence.code=f[[2]]$convergence)))
    # optimized.params contains extra info on the convergence of the gradient fit
    list(chrom=chrom, fittab=fittab, optimized.params=lapply(fits, function(f) f[[2]]))
})

details <- all.data # lapply(all.data, function(d) d$optimized.params)
fits <- do.call(rbind, lapply(all.data, function(d) {
    cbind(chr=d$chrom, d$fittab[which.max(logp)])
}))

save(fits, file=out.fits)
save(details, file=out.details)

if ('snakemake' %in% ls()) {
    sink()
}
