#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
    cat("NOTE: do not specify the final .gz suffix for output files here; it will be created automatically. If '.gz' is included, it will be automatically removed.\n")
    stop("usage: pregenotyping.R integrated_table.tab.gz sc_cigars.tab.gz bulk_cigars.tab.gz scan2_config.yaml out_cigardata.tab.gz out.fdr.rda")
}

int.tab=args[1]
sccigars=args[2]
bulkcigars=args[3]
scan2config=args[4]
out.cigars=sub('.gz$', '', args[5])
out.fdr.rda=args[6]

for (outfile in c(paste0(out.cigars, c('', '.gz', '.tbi')), out.fdr.rda)) {
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
s <- read.integrated.table(s, int.tab, quiet=TRUE)

s <- add.cigar.data(s, sccigars, bulkcigars, quiet=TRUE)
null.sites <- cigar.get.null.sites(s, legacy=TRUE, quiet=TRUE)
null.sites <- null.sites[, .(chr, pos, refnt, altnt, muttype, id.score.x, id.score.y, hs.score.x, hs.score.y)]
colnames(null.sites)[1] <- "#chr"
fwrite(null.sites, file=out.cigars, sep='\t')
out.cigars.gz <- paste0(out.cigars, '.gz')
Rsamtools::bgzip(out.cigars, out.cigars.gz)
Rsamtools::indexTabix(file=out.cigars.gz, format='vcf', comment='#')

s <- compute.fdr.prior.data(s)
fdr.prior.data <- s@fdr.prior.data
save(fdr.prior.data, file=out.fdr.rda)
