#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser()

parser$add_argument('output_table', type='character',
    help='Tab-delimited output table containing all VAF-based calls and signature-rescued called in this batch. This table DOES NOT include additional mutations from --add-muts, if specified.')
parser$add_argument('--sig-homogeneity-test', type='character', metavar='FILE', default=NULL,
    help="Write out results of the signature homogeneity test to FILE. The homogeneity test tests the null hypothesis that high quality somatic mutations from sample X match the batch-wide 'true' mutation signature. One p-value is produced for each sample and mutation type in this batch. The homogeneity test is computed and stored in the SCAN2 output objects in --object regardless of whether this table is written.")
parser$add_argument('--object', nargs=2, action='append', metavar=c('in.rda', 'out.rda'),
    help='This argument requires two values: the first is the input RDA file containing a valid SCAN2 object. The second is an output file (which cannot already exist) to which a new SCAN2 object containing signature-rescued calls will be written. This argument can be specified multiple times to create a batch.')
parser$add_argument('--add-muts', metavar='FILE', type='character', default=NULL,
    help='CSV file containing somatic mutations, one per line, that should be added to the batch for calculating the true mutation specrum. The file must contain at least two columns: muttype (entries in this column can have value "snv" or "indel") and mutsig (values in this column must be SBS96- or ID83-channel values; e.g. ACC:C>A for SBS96 or 3:Del:R:0 for ID83)')
parser$add_argument("--rescue-target-fdr", metavar='FLOAT', type='double', default=0.01,
    help='Target FDR cutoff for calling somatic mutations after rescue adjustments are applied. Similar to SCAN2\'s --target-fdr parameter, this is also not formal FDR control.  Value between 0-1, lower values correspond to more stringent calling.')
parser$add_argument('--threads', default=1, type='integer', metavar='N',
    help='Use N threads to load the SCAN2 objects and perform signature rescue. WARNING: each human-sized SCAN2 object requires ~6 GB RAM per thread.')


# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    # no argument processing
}

args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

out.txt <- args$output_table
files <- args$object
n.cores <- args$threads
add.muts <- args$add_muts
homogeneity.file <- args$sig_homogeneity_test
rescue.target.fdr <- as.numeric(args$rescue_target_fdr)

if (rescue.target.fdr < 0 | rescue.target.fdr > 1)
    stop("--rescue-target-fdr must be a floating point number in [0, 1].")

if (length(files) < 1)
    stop('must supply one or more --object arguments')

in.rdas <- files[,1]
out.rdas <- files[,2]

already.exists <- sapply(c(homogeneity.file, out.txt, out.rdas), file.exists)
if (any(already.exists))
    stop(paste('output file(s) already exist, please delete them first: ',
        c(homogeneity.file, out.txt, out.rdas)[already.exists], collapse='\n'))

suppressMessages(library(scan2))
suppressMessages(library(future))
suppressMessages(library(progressr))
plan(multicore, workers=n.cores)

if (!is.null(add.muts))
    add.muts <- data.table::fread(add.muts)

summary <- scan2::mutsig.rescue(setNames(in.rdas, out.rdas), add.muts=add.muts,
    rescue.target.fdr=rescue.target.fdr)

data.table::fwrite(summary$muts, file=out.txt, sep='\t')

if (!is.null(homogeneity.file)) {
    tests <- do.call(rbind, lapply(1:length(summary$sig.homogeneity.tests), function(i)
        data.frame(sample=names(summary$sig.homogeneity.tests)[i], summary$sig.homogeneity.tests[[i]])))
    data.table::fwrite(tests, file=homogeneity.file, sep='\t')
}

if ('snakemake' %in% ls()) {
    sink()
}
