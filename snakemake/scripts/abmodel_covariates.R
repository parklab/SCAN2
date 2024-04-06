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
            snakemake@input['integrated_table'],
            snakemake@input['abfits'],
            snakemake@input['bedgz'],
            snakemake@params['sc_sample'],
            snakemake@params['genome'],
            snakemake@output['tab'],
            snakemake@output['tabgz'],
            snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 8) {
    stop("usage: spatial_sens_abmodel.R integrated_table.tab.gz abfits.rda regions.bed.gz singlecell_id genome_string out.tab out.tab.gz n_cores")
}

integrated.table.path <- args[1]
abfits.path <- args[2]
bed.path <- args[3]
sc.sample <- args[4]
genome.string <- args[5]
out.tab <- args[6]
out.tab.gz <- args[7]
n.cores <- as.integer(args[8])

for (f in c(out.tab, out.tab.gz))
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))

suppressMessages(library(scan2))
suppressMessages(library(Rsamtools))
suppressMessages(library(future))
suppressMessages(library(progressr))
plan(multicore, workers=n.cores)

ab.fits <- get(load(abfits.path))

sens.regions <- data.table::fread(bed.path)
colnames(sens.regions) <- c('chr', 'start', 'end')
sens.regions$end <- sens.regions$end - 1

gr.sens.regions <- GenomicRanges::GRanges(seqnames=sens.regions$chr,
    ranges=IRanges::IRanges(start=sens.regions$start, end=sens.regions$end),
    seqinfo=genome.string.to.seqinfo.object(genome.string))

# Don't want that many tiles here; each one causes two extra tabix reads.
grs.for.para <- analysis.set.tiling.for.parallelization.helper(
    regions=GenomicRanges::reduce(gr.sens.regions),
    total.tiles=10)

with_progress({
    handlers(handler_newline())
    abmodel <- compute.spatial.sensitivity.abmodel(
        single.cell.id=sc.sample,
        ab.fits=ab.fits,
        integrated.table.path=integrated.table.path,
        grs.for.sens=gr.sens.regions,
        grs.for.parallelization=grs.for.para,
        genome.string=genome.string)
}, enable=TRUE)

cat("Writing results to", out.tab, "\n")
colnames(abmodel)[1] <- "#chr"
data.table::fwrite(abmodel, file=out.tab, sep="\t", col.names=TRUE)

Rsamtools::bgzip(out.tab, out.tab.gz)
Rsamtools::indexTabix(file=out.tab.gz, format='bed', comment='#')

if ('snakemake' %in% ls()) {
    sink()
}
