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
            snakemake@input['config_yaml'],
            snakemake@params['single_cell'],
            snakemake@input['inttab'],
            snakemake@input['abfits'],
            snakemake@output['tab'],
            snakemake@output['tabgz'],
            snakemake@threads,
            snakemake@params['chroms']
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
    stop("usage: compute_ab_ests_and_models.R config.yaml single.cell.ID integrated_table.tab.gz abfits.rda out.tab out.tab.gz n.cores [ chromosomes ]")
}

config.yaml <- args[1]
single.cell <- args[2]
int.tab <- args[3]
abfits <- args[4]
out.tab <- args[5]
out.tab.gz <- args[6]
n.cores <- as.integer(args[7])

chroms <- c()
if (length(args) > 7)
    chroms <- args[-(1:7)]

for (f in c(out.tab, out.tab.gz, paste0(out.tab.gz, '.tbi'))) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

suppressMessages(library(scan2))
suppressMessages(library(future))
suppressMessages(library(progressr))

if (n.cores > 1) {
    plan(multicore, workers=n.cores)
} else {
    plan(sequential)
}

object <- make.scan(config.path=config.yaml, single.cell=single.cell)

# Restrict the analysis set to a set of chromosomes, if specified
if (length(chroms) > 0)
    object@analysis.regions <- object@analysis.regions[seqnames(object@analysis.regions) %in% chroms,]
print(object@analysis.regions)

# Use fewer tiles when run on individual chromosomes to avoid unnecessary
# overhead from futures.
target.tile.number <- 300
if (length(chroms) > 0) {
    target.tile.number <- 100
}

cat('target.tile.number', target.tile.number, '\n')

with_progress({
    handlers(handler_newline())
    object <- run.chunked.pipeline(object=object, int.tab=int.tab, abfits=abfits,
        grs.for.parallelization=analysis.set.tiling.for.parallelization(object, total.tiles=target.tile.number),
        what.to.compute='ab.ests.and.models',
        verbose=FALSE, report.mem=TRUE)
}, enable=TRUE)

print(object@gatk[,.(chr,pos,refnt,altnt,ab,gp.mu,gp.sd,abc.pv,lysis.pv,lysis.beta,mda.pv,mda.beta)])

# hack to make the first line a #-initiated header line for tabix
colnames(object@gatk)[1] <- '#chr'
data.table::fwrite(object@gatk[,.(`#chr`,pos,refnt,altnt,ab,gp.mu,gp.sd,abc.pv,lysis.pv,lysis.beta,mda.pv,mda.beta)],
    file=out.tab, sep='\t', quote=FALSE)
Rsamtools::bgzip(out.tab, out.tab.gz)
Rsamtools::indexTabix(file=out.tab.gz, format='vcf', comment='#')

if ('snakemake' %in% ls()) {
    sink()
}
