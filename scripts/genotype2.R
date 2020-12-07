#!/usr/bin/env Rscript

library(scan2)
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 16)
    stop("usage: Rscript genotype.R mmq60.tab mmq1.tab sc.sample bulk.sample somatic_ab.rda sc_cigars.tab bulk_cigars.tab cigar_tuning.rda outfile.rda fdr fdr_tuning.rda { somatic | spikein } min.sc.alt min.sc.dp min.bulk.dp chromosome")

hmq.file <- args[1]
lmq.file <- args[2]
sc.sample <- args[3]
bulk.sample <- args[4]
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
chromosome <- args[16]

if (file.exists(outfile))
    stop(sprintf("output file %s already exists. please delete it first", outfile))

if (spikein != 'somatic' & spikein != 'spikein')
    stop(sprintf("expected 'somatic' or 'spikein' but got %s", spikein))

gt <- make.scan(single.cell=sc.sample, bulk=bulk.sample)
gt <- read.gatk(gt, path=hmq.file)
gt <- gt[chrom(gt) == chromosome,]
gt <- read.gatk.lowmq(gt, path=lmq.file)
gt <- add.ab.estimates(gt, ab.file)
gt <- compute.models(gt)

sc.cigars <- read.table(sc.cigar.file, stringsAsFactors=FALSE, header=TRUE)
bulk.cigars <- read.table(bulk.cigar.file, stringsAsFactors=FALSE, header=TRUE)
load(cigar.tuning.file)
gt <- add.cigar.data(gt, sc.cigars, bulk.cigars, cigar.training)
gt <- add.static.filters(gt, min.sc.alt=min.sc.alt,
    min.sc.dp=min.sc.dp, min.bulk.dp=min.bulk.dp)

save(gt, file=outfile)
