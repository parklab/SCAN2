#!/usr/bin/env Rscript


library(argparse)

parser <- ArgumentParser()

parser$add_argument("mutations", type="character",
    help='CSV file containing a table of mutations annotated per sample.')
parser$add_argument("germline_control", type='character',
    help='CSV file containing heterozygous SNP sites subjected to leave-1-out SCAN2 calling. Used to estimate somatic sensitivity.')
parser$add_argument("callable_regions", type='character',
    help='Two column tab-delimited text file: a sample column describing each single cell sample and a column with the number of basepairs accessible to SCAN2 analysis.')
parser$add_argument("output_summary", type='character',
    help='Write metadata joined to mutation burden estimates (in CSV format) to this file.')
parser$add_argument('--tag', type='character', metavar='STRING', default=NULL,
    help='Optional identifier tag to add as a first column in both outputs. This can be helpful when merging multiple tables (e.g., from different projects or SNV/indels).')
parser$add_argument("--metadata", type='character', metavar='FILE_PATH',
    help='Optional CSV file containing metadata to which mutation burden estimates will be joined. At minimum, must include a single "sample" column identifying each single cell. This can also specify a subset of samples present in the mutations file.  If no metadata is provided, all samples in the mutations file will be analyzed.')
parser$add_argument("--signature-correction", type='character', metavar='FILE_PATH',
    help='Optional sensitivity correction factors based on mutation signature. Example use includes correcting for the different sensitivity of indel classes in the ID83 indel signature format.')
parser$add_argument("--load-priors", default=FALSE, action='store_true',
    help='[NOT IMPLEMENTED] Read upper bounds on mutation burden as determined by the PRE-GENOTYPING step.')
parser$add_argument("--gbp-per-genome", default=5.845001134, type='double',
    help='Number of haploid basepairs (in billions) per genome. The default value of 5.845001134 corresponds to AUTOSOMES ONLY as determined by GRCh37.')
parser$add_argument("--min-bulk-dp", default=11, type='integer',
    help="Minimum required bulk depth for SCAN2 calling. Must match the values used for calling mutations and generating callable regions.")
parser$add_argument("--min-sc-dp", default=6, type='integer',
    help="Minimum required single cell depth for SCAN2 calling. Must match the values used for calling mutations and generating callable regions.")

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(snakemake@params)
    cat('Intercepted command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

if (file.exists(args$output_summary))
    stop(paste('output file', args$output_summary, 'already exists, please delete it first'))


suppressMessages(library(data.table))  # read large tables like the germline control sites
suppressMessages(library(lme4))
suppressMessages(library(lmerTest)) # adds p-value estimate to lme4 output


cat("Reading mutation table..\n")
muts <- fread(args$mutations)
setkey(muts, sample)
print(muts)

cat("Reading germline control table..\n")
germline <- fread(args$germline_control)
setkey(germline, sample)
print(germline)

cat("Reading callable regions..\n")
callable.regions <- fread(args$callable_regions)
setkey(callable.regions, sample)
print(callable.regions)

if (is.null(args$metadata)) {
    meta <- data.table(sample=unique(muts[['sample']]))
} else {
    cat("Reading metadata..\n")
    meta <- fread(args$metadata)
}
setkey(meta, sample)
print(meta)

if (is.null(args$signature_correction)) {
    muttypes <- unique(c(muts$muttype, germline$muttype))
    correction.factor <- rep(1, length(muttypes))
    names(correction.factor) <- muttypes
} else {
    stop('Correction by signature is not yet implemented.\n')
    cat("Reading signature-based correction factors..\n")
    tmp <- fread(args$signature_correction)
    correction.factor <- tmp[['factor']]
    names(correction.factor) <- tmp[['muttype']]
}
print(correction.factor)




# load the pre-genotyping burden estimate. This should bear on
# sensitivity, since it affects N_T and N_A estimates.
# XXX: Not incorporated for now. The issue is the need to access the
# full SCAN2 output directories.
if (args$load_priors) {
    stop('Importing prior mutation burdens is not yet implemented\n')
    # XXX: donor is no longer a required argument in the metadata table
    donors <- unique(meta[samples,]$donor)
    fcs <- do.call(c, lapply(donors, function(dn) {
        this.samples <- meta$sample[meta$donor==dn & meta$sample %in% samples]
        ret <- lapply(this.samples, function(sn) {
            f <- sprintf("~/ndata1/pta/%s/scansnv_fdr01_noX/indel/%s/fdr_tuning.rda", dn, sn)
            if (!file.exists(f))
                return(list(burden=c(NA,NA)))
            load(f)
            fdr.tuning
        })
        names(ret) <- this.samples
        ret
    }))
    pre.geno.burdens <- sapply(fcs, function(fc) fc$burden[2])
}

samples <- meta[['sample']]

# Compute sensitivity in several ways:
#   1. completely raw sensitivity: fraction of all germline variants that pass
#   2. fraction of germline variants that pass and meet the minimum requirements for being analyzed
#   3. mutation signature-based correction of both (1) and (2). Primarily intended for indels.

# XXX: Start of a more general correction by signature class
correct <- function(muts, correction.factors) {
    mt <- table(muts$muttype[muts$pass])
    sum(mt) / sum(mt/correction.factors[names(mt)])
}



# Get sensitivity estimates from germline variants
sens <- germline[, .(sens=mean(pass),
                     callable.sens=mean(pass[bulk.dp >= args$min_bulk_dp & dp >= args$min_sc_dp])),
                   by=sample]
setkey(sens, sample)

# Get counts of mutations per sample
muttab <- muts[, .(nsom=sum(passA)), by=sample]
setkey(muttab, sample)

# Mash everything together
muttab <- meta[muttab, on='sample'][sens, on='sample', nomatch=0][callable.regions[,.(sample, callable.bp=cbp)], on='sample', nomatch=0]
# Subset back down to whatever was in meta (other files could contain supersets).
muttab <- muttab[sample %in% meta[['sample']]]

# "Callable" estimates are more appropriate when trying to estimate
# mutations per gbp, since low depth regions are ignored. They are
# less useful for performance assessment, where lack of calls in
# regions with low coverage should count against performance.
# The cutoffs of 10 and 5 correspond to SCAN2 defaults.
muttab$callable.burden <- muttab[['nsom']] / muttab[['callable.sens']]
muttab$rate.per.gb <- muttab[['callable.burden']] / muttab[['callable.bp']] * 1e9 / 2
muttab$genome.burden <- muttab[['rate.per.gb']] * args$gbp_per_genome
muttab$genome.sens <- muttab[['nsom']] / muttab[['genome.burden']]

# XXX: --------
# Resurrect at some point to allow signature-based correction
#muttab$burden.id83 <- muttab$nsom / muttab$sens.id83
#muttab$rate.per.gb.id83 <- muttab$burden.id83 / muttab$callable.bp * 1e9 / 2
#muttab$genome.burden.id83 <- round(gbp.to.genome.factor * muttab$rate.per.gb.id83)

# Add optional tag
if (!is.null(args$tag)) {
    muttab <- cbind(args$tag, muttab)
}
fwrite(muttab, file=args$output_summary)
