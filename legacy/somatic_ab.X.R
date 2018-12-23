#!/usr/bin/env Rscript

require('scansnv')
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4)
    stop("usage: somatic_ab.R somatic_positions.txt training.rda fit.rda output.rda\n")

infile <- args[1]
rda <- args[2]
fitfile <- args[3]
outfile <- args[4]

if (file.exists(outfile))
    stop(sprintf("output file %s already exists, please delete it first", outfile))

som.tab <- read.table(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE)
str(som.tab)

load(rda)      # loads 'data'
load(fitfile)  # loads 'fits'
print(fits)

f <- function(chrom) {
    hsnps <- data[data$chr==chrom,]
    fit <- as.data.frame(fits[fits[,1] == chrom, -c(1,6), drop=FALSE])

    pos <- som.tab[som.tab$chr == chrom,]$pos
    if (length(pos) == 0) {
        cat(sprintf("no sites to infer on chr%s\n", chrom))
        return(c())
    }

    cat(sprintf("inferring AB for %d sites on chr%s:%d-%d\n", 
        length(pos), chrom, min(pos), max(pos)))
str(hsnps)
str(fits)
str(fit)
    system.time(z <- infer.gp(ssnvs=data.frame(pos=pos),
                              fit=fit, hsnps=hsnps, chunk=1, flank=2e5))
    cbind(chr=chrom, pos=pos, z)
}

somatic.ab.bychrom <- f("X")

#somatic.ab <- do.call(rbind, somatic.ab.bychrom)
somatic.ab <- somatic.ab.bychrom

save(somatic.ab, file=outfile)
