#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser(description='Prepare somatic mutations for enrichment analysis by filtering recurrent mutations and clustered mutations in the same sample. IMPORTANT: SNVs and indels are filtered separately; so if a SNV and indel are near each other in the same sample, they will not be detected by the cluster filter.')

parser$add_argument('output', type='character',
    help='Output table containing ALL input mutations with a column for each filter and a column named `final.filter` that indicates whether the mutation should be retained for further analysis. SNVs and indels, pass and mutsig rescue mutations are in this table.')
parser$add_argument('--metadata', metavar='FILE', default=NULL,
    help='CSV file mapping single cell sample IDs (as in the sample column of --muts files) to individual IDs (i.e., brain donors). If this file is not specified, then each cell will be treated as though it comes from a different donor. This affects the recurrence filter, which removes any mutation occurring in more than one individual. If a mutation occurs in multiple cells from the same individual, such as a lineage marker, then one of them is retained but this can only be determined if sample->subject metadata is given.')
parser$add_argument('--individual-column', metavar='STRING', default='donor',
    help='Name of the column in --metadata that contains the individual ID.')
parser$add_argument('--sample-column', metavar='STRING', default='sample',
    help='Name of the column in --metadata that contains the sample ID.')
parser$add_argument('--muts', action='append', metavar='FILE', required=TRUE,
    help='CSV file of somatic mutations with at least the following columns: sample, chr, pos, refnt, altnt, muttype, mutsig, pass, rescue. This argument can be specified multiple times to combine tables from multiple runs. Must be specified at least once.')
parser$add_argument("--cluster-filter-bp", metavar='INT', type='integer', default=50,
    help='Remove mutations that are within INT base pairs of the nearest mutation of the same type (i.e., snv or indel) AND in the same sample.')
parser$add_argument("--quiet", action='store_true', default=FALSE,
    help='Do not print any output messages. Useful when output=stdout to pipe into downstream scripts.')

args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

if (length(args$muts) < 1)
    stop('at least one --muts argument must be given')

muttypes <- c('snv', 'indel', 'dbs', 'mnv')

main.output.table <- args$output
quiet.mode <- args$quiet

#print(file_test(op='-f', main.output.table))
#if (file.exists(main.output.table) & file_test(op='-f', main.output.table))
# file_test(-f) DOES NOT work like bash shell -f. Returns TRUE for /dev/stdout.
if (file.exists(main.output.table) & main.output.table != '/dev/stdout' & main.output.table != '/dev/null') # Just a hack
    stop(paste('the output file', main.output.table, 'already exists, please delete it first'))


suppressMessages(library(scan2))

muts <- do.call(rbind, lapply(args$muts, function(mutfile)
    data.table::fread(mutfile)[, .(sample, chr, pos, refnt, altnt, muttype, mutsig, pass, rescue)]))

# Default (no metadata): each cell is treated as coming from a unique subject
sample.to.subject.map <- setNames(unique(muts$sample), unique(muts$sample))
if (!is.null(args$metadata)) {
    meta <- fread(args$metadata)#[sample %in% sample.to.subject.map]
    sample.to.subject.map <- setNames(meta[[args$individual_column]], meta[[args$sample_column]])
    if (!quiet.mode) {
        cat("got metadata sample to subject map:\n")
        print(sample.to.subject.map)
    }
}


muts <- muts[order(chr, pos),]
muts[, subject := sample.to.subject.map[sample]]
muts[, id := paste(chr, pos, refnt, altnt)]

# Filtering for things like recurrence/clustering is done separately
# for SNVs and indels.  Perhaps clustering would be better if the two
# were combined?
for (mt in muttypes) {
    m <- muts[muttype == mt]  # never write to m, just use this to prevent re-filtering for muttype
    if (nrow(m) == 0)
        next

    if (!quiet.mode) cat('Initial statistics:', sum(m$pass), 'passed', sum(m$rescue), 'rescued\n')

    # Remove mutations called in 2 different subjects.
    # This does NOT filter SNVs called more than once in the same
    # individual (=likely lineage marker), in which case one of the
    # multply-called mutations is retained. That filter is applied later.
    if (!quiet.mode) {
        cat('Raw recurrence rates:\n')
        print(table(table(m$id)))
    }
    
    if (!quiet.mode) cat('Recurrence x subject table\n')
    z <- split(m$subject, m$id)
    subjects <- sapply(z, function(v) length(unique(v)))
    recs <- sapply(z, length)
    if (!quiet.mode) print(addmargins(table(recs, subjects)))
    # writing filter field: use the real data.table
    muts[muttype == mt, rec.filter := subjects[id] > 1]


    # nearby mutations in a single sample are more likely to
    # be artifacts. Remove the whole cluster, because it is often
    # true that the entire cluster is caused by the same few reads
    # that probably align poorly or are clipped.
    compute.nearest(muts, by.sample=TRUE)
    muts[, cluster.filter := nearest < args$cluster_filter_bp]
    if (!quiet.mode) {
        cat(sprintf('Flagged %d sites within %d bp in the same sample for removal\n',
            sum(muts$cluster.filter), args$cluster_filter_bp))
        print(addmargins(table(muts$sample, muts$cluster.filter)))
    }

    # if there are any recurrent mutations remaining,
    # then they recur in the same subject and are likely
    # lineage-related true mutations. it's important not to count them as
    # multiple independent mutations or high clonality lineage markers
    # could drive enrichment signals.
    muts[muttype == mt, lineage.filter := duplicated(id)]
    muts[muttype == mt, final.filter := rec.filter | cluster.filter | lineage.filter]
}

data.table::fwrite(muts, file=main.output.table)
