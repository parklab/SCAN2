#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4)
    stop("usage: permute.r n muts.rda callable.bed output.rda")

nperms <- as.integer(args[1])
mutrda <- args[2]
callable.bed <- args[3]
outrda <- args[4]

if (file.exists(outrda))
    stop(sprintf("output file %s already exists, please delete it first", outrda))


library(annotatr)
library(regioneR)

load(mutrda) # loads "somatic"

somatic <- somatic[!is.na(somatic$pass) & somatic$pass,]
cat("making mut granges\n")
str(somatic)
if (nrow(somatic) > 0) {
    muts <- GRanges(
        seqnames=paste0("chr", somatic$chr),
        ranges=IRanges(start=somatic$pos-1, end=somatic$pos))
} else {
    muts <- GRanges()
}

cat("reading callable regions\n")
cf <- read.table(callable.bed, header=F, stringsAsFactors=F)
callable <- GRanges(
    seqnames=paste0('chr', cf[,1]),
    ranges=IRanges(start=cf[,2], cf[,3]))

perms <- lapply(1:nperms, function(i) {
    cat('.')
    randomizeRegions(muts, mask=gaps(callable), allow.overlaps=T, per.chromosome=T)
})
cat('\n')
save(perms, callable, file=outrda)

