#!/usr/bin/env Rscript

require('scansnv')

args <- commandArgs(TRUE)
if (length(args) != 3 & length(args) != 4) {
    stop("Rscript combine.R grid.dir gridn ngrids [output.rda]")
}

grid.dir <- args[1]
gridn <- as.integer(args[2])
ngrids <- as.integer(args[3])

write.to <- NA
if (length(args) == 4)
    write.to <- args[4]

x <- read.fit.one(gridn=gridn, dir=grid.dir, n=ngrids)

if (is.na(write.to)) {
    # the chromosome is ignored here anyway
    b <- build.bounds(list(x), chrs=0)

    # output is: one new boundary per line
    # a pair of lines defines the upper and lower bounds for a parameter
    cat(sprintf("%0.10f\n", b[2]))
    cat(sprintf("%0.10f\n", b[3]))
    cat(sprintf("%0.10f\n", b[4]))
    cat(sprintf("%0.10f\n", b[5]))
    cat(sprintf("%0.10f\n", b[6]))
    cat(sprintf("%0.10f\n", b[7]))
    cat(sprintf("%0.10f\n", b[8]))
    cat(sprintf("%0.10f\n", b[9]))
} else {
    fits <- best.fit(list(x))
    save(fits, file=write.to)
}
