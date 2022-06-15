#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
    cat("NOTE: do not specify the final .gz suffix for output files here; it will be created automatically. If '.gz' is included, it will be automatically removed.\n")
    stop("usage: pregenotyping.R sc_sample_id bulk_sample_id integrated_table.tab.gz sc_cigars.tab.gz bulk_cigars.tab.gz genome out_cigardata.tab.gz")
}

sc.id=args[1]
bulk.id=args[2]
int.tab=args[3]
sccigars=args[4]
bulkcigars=args[5]
genome=args[6]
out.cigars=sub('.gz$', '', args[7])

for (outfile in c(paste0(out.cigars, c('', '.gz', '.tbi')))) {
    if (file.exists(outfile))
        stop(paste0('output file ', outfile, ' already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(Rsamtools))

print(bulk.id)
print(sc.id)
print(genome)

s <- make.scan(sc.id, bulk.id, genome)
s <- read.integrated.table(s, int.tab, quiet=TRUE)
s <- add.cigar.data(s, sccigars, bulkcigars, quiet=TRUE)
null.sites <- cigar.get.null.sites(s, legacy=TRUE, quiet=TRUE)
#null.sites <- null.sites[, .(chr, pos, refnt, altnt, muttype, id.score.x, id.score.y, hs.score.x, hs.score.y)]
null.sites <- null.sites[, .(chr, pos, refnt, altnt, muttype, M.cigars, ID.cigars, HS.cigars, other.cigars, dp.cigars, M.cigars.bulk, ID.cigars.bulk, HS.cigars.bulk, other.cigars.bulk, dp.cigars.bulk)]
colnames(null.sites)[1] <- "#chr"
data.table::fwrite(null.sites, file=out.cigars, sep='\t')
out.cigars.gz <- paste0(out.cigars, '.gz')
Rsamtools::bgzip(out.cigars, out.cigars.gz)
Rsamtools::indexTabix(file=out.cigars.gz, format='vcf', comment='#')
