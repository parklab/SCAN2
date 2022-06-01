#!/usr/bin/env Rscript

if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) {
        ret <- unlist(c(
            snakemake@input['hsnps'],
            snakemake@params['chrom'],
            snakemake@output['rda'],
            snakemake@params['n_tiles'],
            n_cores=snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 10) {
    cat("NOTE: do not specify the final .gz suffix for output files here; it will be created automatically. If '.gz' is included, it will be automatically removed.\n")
    stop("usage: process_gatk.R single.cell.id bulk.id mmq60.tab mmq1.tab phased_all.vcf out.tab")
}

sc.sample=args[1]
bulk.sample=args[2]
gatk=args[3]
gatklowmq=args[4]
phased.vcf=args[5]
out.tab=sub('.gz$', '', args[6])

for (outfile in c(paste0(out.tab, c('', '.gz', '.tbi')))) {
    if (file.exists(outfile))
        stop(paste0('output file ', outfile, ' already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(Rsamtools))

# Maybe one day we'll ask for more config variables (like sc.min.alt, etc.)
# XXX: however, that could come from snakemake@config
#suppressMessages(library(yaml))
#y <- yaml::read_yaml(scan2config)

x <- make.gatk.table(gatk, gatklowmq, phased.vcf, sc.sample, bulk.sample)
colnames(x)[1] <- paste0('#', colnames(x)[1])  # hackish way to begin the header line with #
data.table::fwrite(x, sep='\t', file=out.tab)
out.tab.gz <- paste0(out.tab, '.gz')
Rsamtools::bgzip(file=out.tab, dest=out.tab.gz, overwrite=FALSE)
Rsamtools::indexTabix(file=out.tab.gz, format='vcf', comment='#')

if ('snakemake' %in% ls()) {
    sink()
}
