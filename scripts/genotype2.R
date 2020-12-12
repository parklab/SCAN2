#!/usr/bin/env Rscript

library(scan2)
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 10)
    stop("usage: Rscript genotype.R mmq60.tab mmq1.tab sc.sample bulk.sample somatic_ab.rda outfile.rda min.sc.alt min.sc.dp min.bulk.dp chromosome")

hmq.file <- args[1]
lmq.file <- args[2]
sc.sample <- args[3]
bulk.sample <- args[4]
ab.file <- args[5]
outfile <- args[6]
min.sc.alt <- as.integer(args[7])
min.sc.dp <- as.integer(args[8])
min.bulk.dp <- as.integer(args[9])
chromosome <- args[10]

if (file.exists(outfile))
    stop(sprintf("output file %s already exists. please delete it first", outfile))

gt <- make.scan(single.cell=sc.sample, bulk=bulk.sample)
gt <- read.gatk(gt, path=hmq.file)
gt <- gt[chrom(gt) == chromosome,]
gt <- read.gatk.lowmq(gt, path=lmq.file)
gt <- add.ab.estimates(gt, ab.file)
gt <- compute.models(gt)

#sc.cigars <- read.table(sc.cigar.file, stringsAsFactors=FALSE, header=TRUE)
#bulk.cigars <- read.table(bulk.cigar.file, stringsAsFactors=FALSE, header=TRUE)
#load(cigar.tuning.file)
#gt <- add.cigar.data(gt, sc.cigars, bulk.cigars, cigar.training)
#gt <- add.static.filters(gt, min.sc.alt=min.sc.alt,
    #min.sc.dp=min.sc.dp, min.bulk.dp=min.bulk.dp)

save(gt, file=outfile)
