#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 10) {
    cat("NOTE: do not specify the final .gz suffix for output files here; it will be created automatically\n")
    stop("usage: preprocess.R gatk.tab.gz gatk_lowmq.tab.gz hsnps.tab.gz sc_cigars.tab.gz bulk_cigars.tab.gz scan2_config.yaml out_hsnp_resampling.tab.gz out_hsnp_resampling_details.rda out_cigars.tab.gz out.fdr.rda")
}

gatk=args[1]
gatklowmq=args[2]
hsnps=args[3]
sccigars=args[4]
bulkcigars=args[5]
scan2config=args[6]
out.hsnps=args[7]
out.hsnps.rda=args[8]
out.cigars=args[9]
out.fdr.rda=args[10]

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
s <- read.gatk.lowmq(s, gatklowmq, quiet=TRUE)
s <- add.training.data(s, hsnps, quiet=TRUE)
s <- resample.training.data(s)
resample.data <- s@resampled.training.data

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

s <- compute.fdr.priors(s)
fdr.prior.data <- s@fdr.prior.data
save(fdr.prior.data, file=out.fdr.rda)
