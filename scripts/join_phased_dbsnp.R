#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) {
        ret <- unlist(c(
            snakemake@input['bulk_called_vcf'],
            snakemake@input['phased_hsnps_vcf'],
            snakemake@params['genome'],
            snakemake@output['out_vcf'],
            snakemake@output['out_vcf_gz'],
            snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (!any(length(args) == 5:6)) {
    stop("usage: join_phased_hsnps.R bulk_called.vcf.gz phased_hsnps.vcf.gz genome_string out.vcf out.vcf.gz [n.cores]")
}

bulk.called.vcf <- args[1]
hsnps.vcf <- args[2]
genome <- args[3]
out.vcf <- args[4]
out.vcf.gz <- args[5]

n.cores <- 1
if (length(args) == 6)
    n.cores <- as.integer(args[6])

for (f in c(out.vcf, out.vcf.gz)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(future))
suppressMessages(library(Rsamtools))
plan(multicore, workers=n.cores)

# Currently the chunking used here isn't configurable by user.
results <- join.phased.hsnps(bulk.called.vcf=bulk.called.vcf, hsnps.vcf=hsnps.vcf, genome=genome)

header <- system(paste0("tabix -H ", bulk.called.vcf), intern=TRUE)
# The final line of the header is the #CHROM...
# There will be a genotype column for each sample, but we just want a single
# 'phasedgt' output column.
header[length(header)] <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tphasedgt"

writeLines(text=header, con=out.vcf)
data.table::fwrite(results, file=out.vcf, sep="\t", na="NA", quote=FALSE, append=TRUE)
Rsamtools::bgzip(out.vcf, out.vcf.gz)
Rsamtools::indexTabix(file=out.vcf.gz, format="vcf", comment="#")

if ('snakemake' %in% ls()) {
    sink()
}
