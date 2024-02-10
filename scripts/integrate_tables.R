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
            snakemake@params['config_yaml'],
            snakemake@params['bulk_sample'],
            snakemake@params['snv_min_bulk_dp'],
            snakemake@params['snv_max_bulk_alt'],
            snakemake@params['snv_max_bulk_af'],
            snakemake@params['snv_max_bulk_binom_prob'],
            snakemake@params['indel_min_bulk_dp'],
            snakemake@params['indel_max_bulk_alt'],
            snakemake@params['indel_max_bulk_af'],
            snakemake@params['indel_max_bulk_binom_prob'],
            snakemake@params['genome'],
            snakemake@params['legacy'],
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
if (length(args) < 19) {
    stop("usage: integrate_tables.R mmq60.tab.gz mmq1.tab.gz phased_hets.vcf.gz cross_sample_panel.tab.gz config.yaml bulk_sample_id snv_min_bulk_dp snv_max_bulk_alt snv_max_bulk_af snv_max_bulk_binom_prob indel_min_bulk_dp indel_max_bulk_alt indel_max_bulk_af indel_max_bulk_binom_prob genome_string legacy out.tab out.tab.gz out_resampling_details.rda [n.cores]")
}

mmq60 <- args[1]
mmq1 <- args[2]
phasing <- args[3]
panel <- args[4]
config.yaml <- args[5]
bulk.sample <- args[6]
snv.min.bulk.dp <- as.integer(args[7])
snv.max.bulk.alt <- as.integer(args[8])
snv.max.bulk.af <- as.numeric(args[9])
snv.max.bulk.binom.prob <- as.numeric(args[10])
indel.min.bulk.dp <- as.integer(args[11])
indel.max.bulk.alt <- as.integer(args[12])
indel.max.bulk.af <- as.numeric(args[13])
indel.max.bulk.binom.prob <- as.numeric(args[14])
genome <- args[15]
legacy <- as.logical(args[16])
out.tab <- args[17]
out.tab.gz <- args[18]
out.rda <- args[19]

n.cores <- 1
if (length(args) == 19)
    n.cores <- as.integer(args[20])

for (f in c(out.tab, out.tab.gz, paste0(out.tab.gz, '.tbi'), out.rda)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(yaml))
suppressMessages(library(future))
plan(multicore, workers=n.cores)

if (is.na(panel))
    panel <- NULL

y <- yaml::read_yaml(config.yaml)
sc.samples <- names(read_yaml(config.yaml)$sc_bams)

# Currently the chunking used here isn't configurable by user.
results <- make.integrated.table(mmq60, mmq1, phasing, panel=panel,
    bulk.sample=bulk.sample, sc.samples=sc.samples,
    snv.min.bulk.dp=snv.min.bulk.dp,
    indel.min.bulk.dp=indel.min.bulk.dp,
    snv.max.bulk.alt=snv.max.bulk.alt,
    snv.max.bulk.af=snv.max.bulk.af,
    snv.max.bulk.binom.prob=snv.max.bulk.binom.prob,
    indel.max.bulk.alt=indel.max.bulk.alt,
    indel.max.bulk.af=indel.max.bulk.af,
    indel.max.bulk.binom.prob=indel.max.bulk.binom.prob,
    genome=genome,
    legacy=legacy)

inttab <- results$gatk
write.integrated.table(inttab=inttab, out.tab=out.tab, out.tab.gz=out.tab.gz)

resampling.details <- results$resampling.details
save(resampling.details, file=out.rda)

if ('snakemake' %in% ls()) {
    sink()
}
