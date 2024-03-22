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
            # if the user did not provide a panel, this entry is NULL and c() drops
            # NULL values silently. replace with NA and handle it later.
            ifelse(is.null(snakemake@params[['cross_sample_panel']]), NA,
                snakemake@params[['cross_sample_panel']]),
            snakemake@input['config_yaml'],
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
if (length(args) < 8) {
    cat('if no cross sample panel was created, use the string "NA" to skip loading\n')
    stop("usage: integrate_tables.R mmq60.tab.gz mmq1.tab.gz phase_info.tab.gz cross_sample_panel.tab.gz config.yaml out.tab out.tab.gz out_resampling_details.rda [n.cores]")
}

print(args)

mmq60 <- args[1]
mmq1 <- args[2]
phasing <- args[3]
panel <- args[4]
config.yaml <- args[5]
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

if (is.na(panel) | panel == "NA")
    panel <- NULL

# This SCAN2 object is not for use, it is just to parse the config file.
# Keep in mind that it makes _no sense_ to have an object here because
# SCAN2 objects represent one single cell while the integrated table
# represents all cells from an individual.
dummy.object <- make.scan(config.path=config.yaml)

# Currently the chunking used here isn't configurable by user.
results <- make.integrated.table(dummy.object, mmq60, mmq1, phasing, panel=panel)

inttab <- results$gatk
write.integrated.table(inttab=inttab, out.tab=out.tab, out.tab.gz=out.tab.gz)

resampling.details <- results$resampling.details
save(resampling.details, file=out.rda)

if ('snakemake' %in% ls()) {
    sink()
}
