library(scan2)

load(snakemake@input[['gtobject']])

sctab <- read.table(snakemake@input[['sc_cigars']],
    header=TRUE, stringsAsFactors=FALSE)
bulktab <- read.table(snakemake@input[['bulk_cigars']],
    header=TRUE, stringsAsFactors=FALSE)

gt <- gt[chrom(gt) == snakemake@wildcards[['chr']],]
gt <- add.cigar.data(gt, sc.cigars=sctab, bulk.cigars=bulktab)

save(gt, file=snakemake@output[[1]])
