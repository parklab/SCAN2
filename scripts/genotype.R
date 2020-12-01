#!/usr/bin/env Rscript

require('scansnv')
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 15)
    stop("usage: Rscript genotype.R mmq60.tab mmq1.tab sc.sample bulk.sample somatic_ab.rda sc_cigars.tab bulk_cigars.tab cigar_tuning.rda outfile.rda fdr fdr_tuning.rda { somatic | spikein } min.sc.alt min.sc.dp min.bulk.dp")

hmq.file <- args[1]
lmq.file <- args[2]
sc.sample <- make.names(args[3])
bulk.sample <- make.names(args[4])
ab.file <- args[5]
sc.cigar.file <- args[6]
bulk.cigar.file <- args[7]
cigar.tuning.file <- args[8]
outfile <- args[9]
fdr <- as.numeric(args[10])
fdr.tuning.file <- args[11]
spikein <- args[12]
min.sc.alt <- as.integer(args[13])
min.sc.dp <- as.integer(args[14])
min.bulk.dp <- as.integer(args[15])

if (file.exists(outfile))
    stop(sprintf("output file %s already exists. please delete it first", outfile))

if (spikein != 'somatic' & spikein != 'spikein')
    stop(sprintf("expected 'somatic' or 'spikein' but got %s", spikein))

load(cigar.tuning.file)  # loads cigar.emp.score(), cigar.training
load(ab.file, verbose=T) # loads ab
somatic.ab <- ab
str(somatic.ab)


# Step 1: read a few rows just to get column names
hmq <- read.table(hmq.file, header=T, stringsAsFactors=F,
    colClasses=c(chr='character'), nrow=10)
tot.cols <- ncol(hmq)
sc.idx <- which(colnames(hmq) == sc.sample)
bulk.idx <- which(colnames(hmq) == bulk.sample)
cat("Using data:\n")
for (i in 1:ncol(hmq)) {
    s <- ifelse(i == sc.idx, "[SC]", ifelse(i == bulk.idx, "[BULK]", ""))
    cat(sprintf("%8s %s\n", s, colnames(hmq)[i]))
}


# Step 2: really read the tables in, but only the relevant columns
#hmq <- read.table(hmq.file, sep='\t', header=T, stringsAsFactors=F)
cols.to.read <- rep("NULL", tot.cols)
# First 7 are chr, pos, dbsnp, refnt, altnt, mq, mqrs
cols.to.read[1:7] <- c('character', 'integer', rep('character', 3), 'numeric', 'numeric')
# Read 3 columns for the single cell, 3 columns for bulk
cols.to.read[sc.idx + 0:2] <- c('character', 'integer', 'integer')
cols.to.read[bulk.idx + 0:2] <- c('character', 'integer', 'integer')
cat("reading mmq60 GATK table..\n")
cat(sprintf("reading %d columns from %s..\n",
    sum(is.na(cols.to.read) | cols.to.read != 'NULL'), hmq.file))
hmq <- read.table(hmq.file, header=T, stringsAsFactors=F,
    colClasses=cols.to.read)
new.sc.idx <- which(colnames(hmq) == sc.sample)
new.bulk.idx <- which(colnames(hmq) == bulk.sample)
str(hmq)

cat("reading mmq1 GATK table..\n")
cat(sprintf("reading %d columns from %s..\n",
    sum(is.na(cols.to.read) | cols.to.read != 'NULL'), lmq.file))
lmq <- read.table(lmq.file, sep='\t', header=T, stringsAsFactors=F,
    colClasses=cols.to.read)
str(lmq)



sc.cigar <- read.table(sc.cigar.file, sep='\t', stringsAsFactors=F, header=T)
bulk.cigar <- read.table(bulk.cigar.file, sep='\t', stringsAsFactors=F, header=T)

load(fdr.tuning.file) # loads 'fdr.tuning'

gt <- genotype.somatic(gatk=hmq, gatk.lowmq=lmq,
    sc.idx=new.sc.idx, bulk.idx=new.bulk.idx,
    sites.with.ab=somatic.ab, fdr.tuning=fdr.tuning,
    sc.cigars=sc.cigar, bulk.cigars=bulk.cigar,
    cigar.training=cigar.training, cigar.emp.score=cigar.emp.score,
    target.fdr=fdr, spikein=spikein == 'spikein',
    min.sc.alt=min.sc.alt, min.bulk.dp=min.bulk.dp, min.sc.dp=min.sc.dp)

#gt$somatic <- get.3mer(gt$somatic)

save(gt, file=outfile)
