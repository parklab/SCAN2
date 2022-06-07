#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 8) {
    cat("NOTE: do not specify the final .gz suffix for output files here; it will be created automatically. If '.gz' is included, it will be automatically removed.\n")
    stop("usage: pregenotyping.R sc_sample_id bulk_sample_id integrated_table.tab.gz sc_cigars.tab.gz bulk_cigars.tab.gz scan2_config.yaml out_cigardata.tab.gz out.fdr.rda")
}

sc.id=args[1]
bulk.id=args[2]
int.tab=args[3]
sccigars=args[4]
bulkcigars=args[5]
scan2config=args[6]
out.cigars=sub('.gz$', '', args[7])
out.fdr.rda=args[8]

for (outfile in c(paste0(out.cigars, c('', '.gz', '.tbi')), out.fdr.rda)) {
    if (file.exists(outfile))
        stop(paste0('output file ', outfile, ' already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(yaml))
suppressMessages(library(Rsamtools))

y <- yaml::read_yaml(scan2config)

print(bulk.id)
print(sc.id)
print(y$genome)

s <- make.scan(sc.id, bulk.id, y$genome)
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
