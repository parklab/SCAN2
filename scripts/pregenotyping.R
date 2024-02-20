#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
    cat("NOTE: do not specify the final .gz suffix for output files here; it will be created automatically. If '.gz' is included, it will be automatically removed.\n")
    stop("usage: pregenotyping.R config.yaml sc_sample_id integrated_table.tab.gz sc_cigars.tab.gz bulk_cigars.tab.gz out_cigardata.tab.gz")
}

config.yaml=args[1]
sc.id=args[2]
int.tab=args[3]
sccigars=args[4]
bulkcigars=args[5]
out.cigars=sub('.gz$', '', args[6])

for (outfile in c(paste0(out.cigars, c('', '.gz', '.tbi')))) {
    if (file.exists(outfile))
        stop(paste0('output file ', outfile, ' already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(Rsamtools))

print(sc.id)

s <- make.scan(config.path=config.yaml, single.cell=sc.id)
s <- read.integrated.table(s, int.tab, quiet=TRUE)
s <- add.cigar.data(s, sccigars, bulkcigars, quiet=TRUE)
null.sites <- cigar.get.null.sites(s, legacy=TRUE, quiet=TRUE)
null.sites <- null.sites[, .(chr, pos, refnt, altnt, muttype, M.cigars, ID.cigars, HS.cigars, other.cigars, dp.cigars, M.cigars.bulk, ID.cigars.bulk, HS.cigars.bulk, other.cigars.bulk, dp.cigars.bulk)]
colnames(null.sites)[1] <- "#chr"
data.table::fwrite(null.sites, file=out.cigars, sep='\t')
out.cigars.gz <- paste0(out.cigars, '.gz')
Rsamtools::bgzip(out.cigars, out.cigars.gz)
Rsamtools::indexTabix(file=out.cigars.gz, format='vcf', comment='#')
