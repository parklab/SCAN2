#!/usr/bin/env Rscript

require('scansnv')
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3)
    stop("usage: get_hsnp_positions.R training_hsnps.rda output.txt")


rda.file <- args[1]
outfile <- args[2]
n <- as.integer(args[3])

if (file.exists(outfile))
    stop(sprintf("output file %s already exists. please delete it first", outfile))

load(rda.file)  # must contain a data frame named "data"

hsnp.positions <- data[sample(nrow(data), size=n),]

write.table(hsnp.positions[,c('chr', 'pos')], file=outfile, sep='\t',
    col.names=TRUE, row.names=FALSE, quote=FALSE)
