#!/usr/bin/env Rscript

require('scansnv')
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 5)
    stop("usage: get_somatic_positions.R gatk_mmq60.table gatk_mmq1.table sc.sample bulk.sample output.txt")

hmq.file <- args[1]
lmq.file <- args[2]
sc.sample <- make.names(args[3])
bulk.sample <- make.names(args[4])
outfile <- args[5]

if (file.exists(outfile))
    stop(sprintf("output file %s already exists. please delete it first", outfile))

hmq <- read.table(hmq.file, sep="\t", stringsAsFactors=F, header=TRUE)
lmq <- read.table(lmq.file, sep="\t", stringsAsFactors=F, header=TRUE)

sc.idx <- which(colnames(hmq) == sc.sample)
bulk.idx <- which(colnames(hmq) == bulk.sample)
cat("Using data:\n")
for (i in 1:ncol(hmq)) {
    s <- ifelse(i == sc.idx, "[SC]", ifelse(i == bulk.idx, "[BULK]", ""))
    cat(sprintf("%8s %s\n", s, colnames(hmq)[i]))
}

som.positions <- genotype.somatic(hmq, lmq,
    sc.idx=sc.idx, bulk.idx=bulk.idx, sites.only=TRUE)$somatic

write.table(som.positions[,c('chr', 'pos')], file=outfile, sep='\t',
    col.names=TRUE, row.names=FALSE, quote=FALSE)
