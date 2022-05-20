#!/bin/bash

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != ) {
    stop("usage: preprocess.R gatk.tab.gz hsnps.tab.gz sc_cigars.tab.gz bulk_cigars.tab.gz scan2_config.yaml out.tab.gz out_details.rda"
}

gatk=args[1]
hsnps=args[2]
sccigars=args[3]
bulkcigars=args[4]
scan2config=args[5]
out.hsnps=args[6]
out.hsnps.rda=args[7]
out.cigars=args[8]
out.fdr.rda=args[9]

for (outfile in c(paste0(out.hsnps, c('', '.gz', '.tbi')), paste0(out.cigars, c('', '.gz', '.tbi')), out.hsnps.rda, out.fdr.rda)) {
    if (file.exists(outfile))
        stop(paste0('output file ', outfile, ' already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(yaml))
suppressMessages(library(Rsamtools))

y <- yaml::read_yaml(scan2config)

bulk_id <- y$bulk_sample
sc_id <- names(y$sc_bams)[1]
print(bulk_id)
print(sc_id)

# This reads GATK, the hSNPs tab, and the SCAN2 config file and uses
# SCAN2 methods to downsample hSNPs. The resulting table is written to
# stdout for bgzip and tabix.
# NOTE: the genome type in make.scan doesn't matter at all
s <- make.scan(sc_id, bulk_id, 'hs37d5')
s <- add.static.filter.params(s, scan2config)
s <- read.gatk(s, gatk, quiet=TRUE, add.mutsig=FALSE)
s <- add.training.data(s, '$hsnps', quiet=TRUE)
s <- resample.training.data(s)

save(resample.data, file=out.hsnps.rda)
out <- s@gatk[training.site==TRUE,.(chr, pos, refnt, altnt, training.hap1, training.hap2, training.phgt, resampled.training.site)]
colnames(out) <- c('#chr', 'pos', 'refnt', 'altnt', 'hap1', 'hap2', 'phgt', 'resampled')
fwrite(out, file=out.hsnps, sep='\t')
Rsamtools::bgzip(out.hsnps, paste0(out.hsnps, '.gz'))
system(paste('tabix -p vcf -S 1', paste0(out.hsnps, '.gz')))

s <- add.cigar.data(s, sccigars, bulkcigars, quiet=TRUE)
null.sites <- cigar.get.null.sites(s, legacy=TRUE, quiet=TRUE)
null.sites <- null.sites[, .(chr, pos, refnt, altnt, muttype, id.score.x, id.score.y, hs.score.x, hs.score.y)]
colnames(null.sites)[1] <- "#chr"
fwrite(null.sites, file=out.cigars, sep='\t')
Rsamtools::bgzip(out.cigars, paste0(out.cigars, '.gz'))
system(paste('tabix -p vcf -S 1', paste0(out.cigars, '.gz')))

