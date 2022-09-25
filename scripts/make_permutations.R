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
            snakemake@input['muts'],
            snakemake@params['sc_sample'],
            snakemake@input['callable'],
            snakemake@input['genome_file'],
            snakemake@params['genome_string'],
            snakemake@params['muttype'],
            snakemake@params['passtype'],
            snakemake@params['n_permutations'],
            snakemake@params['generation_param'],
            snakemake@output['rda'],
            snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
print(args)
if (length(args) < 10) {
    cat('mutations.csv may contain mutations from multiple single cells so long as a `sample` column is provided. It must also include columns: chr, pos, refnt, altnt, pass, rescue.\n')
    cat("N.B.: genome_string must be a supported SCAN2 genome.\n")
    cat("generation_param - for SNVs, the number of random sites to generate per iteration; for indels, a base multiplier that is linearly related to the number of random sites generated. Ideally, each iteration should generate >1 permutation of the input; for reasonable performance, 10s or 100s of permutations should be generated per iteration. However, using very large values of this parameter can cause excessive memory usage.\n")
    stop("usage: make_permutations.R mutations.csv sc_sample callable_regions.bed genome_string bedtools_genomefile.txt {snv|indel} {pass|rescue} n_permutations out.rda generation_param [n.cores]")
}

muts.path <- args[1]
sc.sample <- args[2]
callable.path <- args[3]
genome.string <- args[4]
genomefile.path <- args[5]
muttype <- args[6]
passtype <- args[7]
n.permutations <- as.integer(args[8])
generation.param <- as.numeric(args[9])
out.rda <- args[10]

n.cores <- 1
if (length(args) == 11)
    n.cores <- as.integer(args[11])

if (muttype != 'snv' & muttype != 'indel')
    stop("muttype must be either 'snv' or 'indel' (case-sensitive)")
mt <- muttype

if (passtype != 'pass' & passtype != 'rescue')
    stop("passtype must be either 'pass' or 'rescue' (case-sensitive)")

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))


suppressMessages(library(scan2))
suppressMessages(library(future))
suppressMessages(library(progressr))
plan(multicore, workers=n.cores)

muts <- data.table::fread(muts.path)
muts <- muts[sample == sc.sample & muttype == mt & (pass == TRUE | (passtype == 'rescue' & rescue == TRUE))]
cat('Got', nrow(muts), 'mutations of type', mt, '\n')

progressr::with_progress({
    # handler_newline causes alot of printing, but it's log-friendly
    progressr::handlers(progressr::handler_newline())
    permdata <- make.permuted.mutations(muts=muts,
        sc.sample=sc.sample,
        callable.bed=callable.path, genome=genome.string,
        genome.file=genomefile.path, muttype=muttype,
        n.permutations=n.permutations, n.chunks=ceiling(n.permutations/100),
        quiet=TRUE,
        # CAUTION!! This hack assumes that make.permuted.mutations
        # always considers SNVs and indels separately, meaning that
        # one of these two generation.param settings will be
        # ignored. The two generation.params ARE NOT interpreted
        # the same way and are not interchangeable.
        snv.N=generation.param, indel.K=generation.param)
}, enable=TRUE)

save(permdata, file=out.rda, compress=FALSE)

if ('snakemake' %in% ls()) {
    sink()
}
