#!/usr/bin/env Rscript

require('scansnv')
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 11)
    stop("usage: Rscript genotype.R mmq60.tab mmq1.tab sc.sample bulk.sample somatic_ab.rda somatic_cigars.rda hsnp_cigars.rda outfile.rda fdr fdr_tuning.rda { somatic | spikein }")

hmq.file <- args[1]
lmq.file <- args[2]
sc.sample <- make.names(args[3])
bulk.sample <- make.names(args[4])
ab.file <- args[5]
somatic.cigar.file <- args[6]
hsnp.cigar.file <- args[7]
outfile <- args[8]
fdr <- as.numeric(args[9])
fdr.tuning.file <- args[10]
spikein <- args[11]

if (file.exists(outfile))
    stop(sprintf("output file %s already exists. please delete it first", outfile))

if (spikein != 'somatic' & spikein != 'spikein')
    stop(sprintf("expected 'somatic' or 'spikein' but got %s", spikein))

load(ab.file, verbose=T) # loads ab
somatic.ab <- ab
str(somatic.ab)

cat("reading mmq60 GATK table..\n")
hmq <- read.table(hmq.file, sep='\t', header=T, stringsAsFactors=F)
cat("reading mmq1 GATK table..\n")
lmq <- read.table(lmq.file, sep='\t', header=T, stringsAsFactors=F)
somatic.cigar <- read.table(somatic.cigar.file, sep='\t', stringsAsFactors=F, header=T)
hsnp.cigar <- read.table(hsnp.cigar.file, sep='\t', stringsAsFactors=F, header=T)

sc.idx <- which(colnames(hmq) == sc.sample)
bulk.idx <- which(colnames(hmq) == bulk.sample)
cat("Using data:\n")
for (i in 1:ncol(hmq)) {
    s <- ifelse(i == sc.idx, "[SC]", ifelse(i == bulk.idx, "[BULK]", ""))
    cat(sprintf("%8s %s\n", s, colnames(hmq)[i]))
}

load(fdr.tuning.file) # loads 'fdr.tuning'

gt <- genotype.somatic(gatk=hmq, gatk.lowmq=lmq,
    sc.idx=sc.idx, bulk.idx=bulk.idx,
    sites.with.ab=somatic.ab, fdr.tuning=fdr.tuning,
    somatic.cigars=somatic.cigar, hsnp.cigars=hsnp.cigar,
    target.fdr=fdr, spikein=spikein == 'spikein')

#gt$somatic <- get.3mer(gt$somatic)

save(gt, file=outfile)
