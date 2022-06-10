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
            snakemake@params['cross_sample_panel'],
            snakemake@params['bulk_sample'],
            snakemake@params['genome'],
            snakemake@output['tab'],
            snakemake@output['tabgz'],
            snakemake@output['details'],
            snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 9) {
    stop("usage: integrate_tables.R mmq60.tab.gz mmq1.tab.gz phased_hets.vcf.gz cross_sample_panel.tab.gz bulk_sample_id genome_string out.tab out.tab.gz out_resampling_details.rda [n.cores]")
}

mmq60 <- args[1]
mmq1 <- args[2]
phasing <- args[3]
panel <- args[4]
bulk.sample <- args[5]
genome <- args[6]
out.tab <- args[7]
out.tab.gz <- args[8]
out.rda <- args[9]

n.cores <- 1
if (length(args) == 10)
    n.cores <- as.integer(args[10])

for (f in c(out.tab, out.tab.gz, paste0(out.tab.gz, '.tbi'), out.rda)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(future))
plan(multicore, workers=n.cores)

# Currently the chunking used here isn't configurable by user.
results <- make.integrated.table(mmq60, mmq1, phasing, panel=panel, bulk.sample, genome)

inttab <- results$gatk
write.integrated.table(inttab=inttab, out.tab=out.tab, out.tab.gz=out.tab.gz)

resampling.details <- results$resampling.details
save(resampling.details, file=out.rda)

if ('snakemake' %in% ls()) {
    sink()
}
