#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (!(length(args) %in% 4:5))
    stop("usage: make_panel.R input.tab.gz metadata.csv scan2config.yaml output.tab [n.cores]")

inf <- args[1]
metaf <- args[2]
configf <- args[3]
outf <- args[4]
n.cores <- 1
if (length(args) == 5)
    n.cores <- as.integer(args[5])
outf.gz <- paste0(outf, '.gz')


for (outfile in paste0(outf, c('', '.gz', '.gz.tbi'))) {
    if (file.exists(outfile))
        stop(sprintf("output file %s already exists, please delete it first", outf))
}

library(scan2)
library(Rsamtools)
library(future)
library(future.apply)
library(progressr)

if (n.cores > 1) {
    plan(multicore, workers=n.cores)
} else {
    plan(sequential)
}


# These panels can be extremely large, so need to read in only the
# columns that really matter.
# Read the first 10 rows, just to get the header.
# Format is: chrom, pos, dbsnp, refnt, altnt, nalleles, N*(gt, ref, alt)
# where N* means repeat the 3 column pattern N times (N=sample number).
# Currently, this just skips the genotype and ref count columns for each
# sample.
read.alts.for.samples <- function(path, region, meta) {
    tf <- Rsamtools::TabixFile(path)
    open(tf)
    # header is not parsed, it's a string of column names separated by tabs
    header <- read.tabix.header(tf)
    close(tf)

    cols <- strsplit(header, '\t')[[1]]

    # All columns past the first 6 are integer counts of alt reads
    cols.to.read <- rep("integer", length(cols))
    cols.to.read[1:6] <- c('character', 'integer', 'character', 'character',
        'character', 'integer')

    # only read in samples that are listed in the metadata
    sample.ids <- c()
    for (i in 7:length(cols)) {
        if (!(cols[i] %in% meta$sample)) {
            cols.to.read[i] <- 'NULL'
            sample.ids <- c(sample.ids, cols[i])
        }
    }

    # Only read the alt counts for samples in the metadata table.
    read.tabix.data(path=path, region=region, header=header,
        colClasses=cols.to.read)
}
    
# Counts the number of unique cells and unique donors that support each mutation
# dmap - a named vector that maps sample ID -> donor ID
make.panel <- function(df, dmap, bulks) {
    if (nrow(df) == 0) {
        return(data.table::data.table(df[,1:6], unique.donors=integer(0),
            unique.cells=integer(0), unique.bulks=integer(0),
            max.out=integer(0), sum.out=integer(0), sum.bulk=integer(0)))
    }

    bulk.idxs <- which(colnames(df) %in% bulks)

    meta.idxs <- 1:6
    sc.idxs <- setdiff(1:ncol(df), c(1:6, bulk.idxs))
    unique.cells <- rowSums(df[,..sc.idxs, drop=FALSE] > 0)
    unique.donors <- rowSums(do.call(cbind, lapply(unique(dmap), function(dn) {
        sns <- names(dmap)[dmap==dn]
        idxs <- which(colnames(df) %in% sns)
        rowSums(df[,..idxs,drop=FALSE]) > 0
    })))
    unique.bulks <- rowSums(df[,..bulk.idxs,drop=FALSE] > 0)

    # remove the entry with the maximum support
    outs <- apply(df[,-..meta.idxs,drop=FALSE], 1, function(row) {
        r <- row[-which.max(row)]
        c(max(r), sum(r))
    })
    max.out <- outs[1,]
    sum.out <- outs[2,]
    sum.bulk <- rowSums(df[,..bulk.idxs,drop=FALSE])
    ret <- data.table(df[,..meta.idxs])
    ret[, c('unique.donors', 'unique.cells', 'unique.bulks', 'max.out', 'sum.out', 'sum.bulk') :=
        list(unique.donors, unique.cells, unique.bulks, max.out, sum.out, sum.bulk)]
    ret
}


# First load metadata. We will not read count data for any samples not
# in the metadata table.
cat(sprintf("Loading sample metadata %s..\n", metaf))
meta <- data.table::fread(metaf)
cat("Bulks: ")
print(meta[amp=='bulk']$sample)

# dummy.object only used to build GRanges intervals for chunked pipeline
config <- scan2::read.config(configf)
config$bulk_sample <- 'PLACEHOLDER'
config$sex <- 'male'
dummy.object <- scan2::make.scan(config=config)
grs <- scan2::analysis.set.tiling.for.parallelization(object=dummy.object)

# map donor <-> sample IDs
dmap <- meta$donor   
names(dmap) <- meta$sample

cat('Starting integrated table pipeline on', length(grs), 'chunks.\n')
cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

progressr::with_progress({
    handlers(handler_newline())
    p <- progressr::progressor(along=1:length(grs))
    p(amount=0, class='sticky', scan2::perfcheck(print.header=TRUE))
    xs <- future.apply::future_lapply(1:length(grs), function(i) {
        gr <- grs[i,]

        pc <- scan2::perfcheck(paste('read.alts.for.samples', i),
            tb <- read.alts.for.samples(path=inf, region=gr, meta=meta))
        p(class='sticky', amount=0, pc)

        pc <- scan2::perfcheck(paste('make.panel', i),
            ret <- make.panel(tb, dmap, meta[amp=='bulk']$sample))
        p(class='sticky', amount=0, pc)

        p()
        ret
    })
}, enable=TRUE)

pon <- data.table::rbindlist(xs)

cat(sprintf("Saving PON to %s..\n", outf))
colnames(pon)[1] <- paste0('#', colnames(pon)[1])
data.table::fwrite(pon, file=outf, sep='\t', quote=FALSE)
Rsamtools::bgzip(outf, outf.gz)
Rsamtools::indexTabix(file=outf.gz, format='vcf', comment='#')
