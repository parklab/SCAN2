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
            snakemake@input['mmq60'],
            snakemake@input['mmq1'],
            snakemake@input['phasing'],
            snakemake@params['bulk_sample'],
            snakemake@params['genome'],
            snakemake@output['tab'],
            snakemake@output['tabgz'],
            snakemake@output['details']
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 8) {
    stop("usage: integrate_tables.R mmq60.tab.gz mmq1.tab.gz phased_hets.vcf.gz bulk_sample_id genome_string out.tab out.tab.gz out_resampling_details.rda [n.cores]")
}

mmq60 <- args[1]
mmq1 <- args[2]
phasing <- args[3]
bulk.sample <- args[4]
genome <- args[5]
out.tab <- args[6]
out.tab.gz <- args[7]
out.rda <- args[8]

n.cores <- 1
if (length(args) == 9)
    n.cores <- as.integer(args[9])

for (f in c(out.tab, out.tab.gz, paste0(out.tab.gz, '.tbi'), out.rda)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(future))
plan(multicore, workers=n.cores)

# Currently the chunking used here isn't configurable by user.
results <- make.integrated.table(mmq60, mmq1, phasing, bulk.sample, genome)

inttab <- results$gatk
colnames(inttab)[1] <- paste0('#', colnames(inttab)[1]) # hack to comment out the header
data.table::fwrite(inttab, file=out.tab, sep='\t', na='NA', quote=FALSE)
Rsamtools::bgzip(out.tab, out.tab.gz)
Rsamtools::indexTabix(file=out.tab.gz, format='vcf', comment='#')

resampling.details <- results$resampling.details
save(resampling.details, file=out.rda)

if ('snakemake' %in% ls()) {
    sink()
}
