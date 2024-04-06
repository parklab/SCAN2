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
            snakemake@input['sccigars'],
            snakemake@input['bulkcigars'],
            snakemake@input['trainingdata'],
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
if (length(args) < 9) {
    stop("usage: compute_ab_ests_and_models.R config.yaml single.cell.ID integrated_table.tab.gz sccigars.tab.gz bulkcigars.tab.gz cigardata.tab.gz out.tab out.tab.gz n.cores [ chromosomes ]")
}

config.yaml <- args[1]
single.cell <- args[2]
int.tab <- args[3]
sccigars.path <- args[4]
bulkcigars.path <- args[5]
trainingcigars.path <- args[6]
out.tab <- args[7]
out.tab.gz <- args[8]
n.cores <- as.integer(args[9])

chroms <- c()
if (length(args) > 9)
    chroms <- args[-(1:9)]

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
    object <- run.chunked.pipeline(object=object, int.tab=int.tab,
        sccigars=sccigars.path, bulkcigars=bulkcigars.path, trainingcigars=trainingcigars.path,
        grs.for.parallelization=analysis.set.tiling.for.parallelization(object, total.tiles=target.tile.number),
        what.to.compute='excess.cigar',
        verbose=FALSE, report.mem=TRUE)
}, enable=TRUE)

print('hi')
print(object@gatk)

# hack to make the first line a #-initiated header line for tabix
colnames(object@gatk)[1] <- '#chr'
data.table::fwrite(object@gatk[,.(`#chr`,pos,refnt,altnt,muttype,id.score,hs.score)],
    file=out.tab, sep='\t', quote=FALSE)
Rsamtools::bgzip(out.tab, out.tab.gz)
Rsamtools::indexTabix(file=out.tab.gz, format='vcf', comment='#')

if ('snakemake' %in% ls()) {
    sink()
}
