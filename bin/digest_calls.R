#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser(description='Prepare somatic mutations for enrichment analysis by filtering recurrent mutations and clustered mutations in the same sample. IMPORTANT: SNVs and indels are filtered separately; so if a SNV and indel are near each other in the same sample, they will not be detected by the cluster filter.')

parser$add_argument('output', type='character',
    help='Output table containing only final, filtered mutations. SNVs and indels, pass and rescue mutations are in this table.')
parser$add_argument('--metadata', metavar='FILE', default=NULL,
    help='CSV file mapping single cell sample IDs (as in the sample column of --muts files) to individual IDs (i.e., brain donors). If this file is not specified, then each cell will be treated as though it comes from a different donor. This affects the recurrence filter, which removes any mutation occurring in more than one individual. If a mutation occurs in multiple cells from the same individual, such as a lineage marker, then one of them is retained but this can only be determined if sample->subject metadata is given.')
parser$add_argument('--muts', action='append', metavar='FILE', required=TRUE,
    help='CSV file of somatic mutations with at least the following columns: sample, chr, pos, refnt, altnt, muttype, mutsig, pass, rescue. This argument can be specified multiple times to combine tables from multiple runs. Must be specified at least once.')
parser$add_argument("--cluster-filter-bp", metavar='INT', type='integer', default=50,
    help='Remove mutations that are within INT base pairs of the nearest mutation of the same type (i.e., snv or indel) AND in the same sample.')
parser$add_argument("--separate-files", type='character', metavar='PREFIX', default=NULL,
    help='Create separate output tables for each combination of pass vs. rescue and snv vs. indel. These tables will be named PREFIX_{pass|pass_and_rescue}_{snv|indel}.txt')

args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

if (length(args$muts) < 1)
    stop('at least one --muts argument must be given')

muttypes <- c('snv', 'indel')

main.output.table <- args$output

output.tables <- paste0(paste(args$separate_files, rep(muttypes, 2), rep(c('pass', 'pass_and_rescue'), each=2), sep='_'), '.txt')
names(output.tables) <- paste(rep(muttypes, 2), rep(c('pass', 'pass_and_rescue'), each=2), sep='_')

check.tables <- main.output.table
if (!is.null(args$separate.files))
    check.tables <- c(check.tables, output.tables)

already.exists <- sapply(check.tables, file.exists)
if (any(already.exists))
    stop(paste('the following output files already exist, please delete them first:', check.tables[already.exists], collapse='\n'))


suppressMessages(library(scan2))

all.muts <- do.call(rbind, lapply(args$muts, function(mutfile)
    data.table::fread(mutfile)[, .(sample, chr, pos, refnt, altnt, muttype, mutsig, pass, rescue)]))

# Default (no metadata): each cell is treated as coming from a unique subject
sample.to.subject.map <- setNames(unique(all.muts$sample), unique(all.muts$sample))
if (!is.null(args$metadata)) {
    meta <- fread(args$metadata)
    sample.to.subject.map <- setNames(meta$subject, meta$sample)
    cat("got metadata sample to subject map:\n")
    print(sample.to.subject.map)
}


all.muts <- all.muts[order(chr, pos),]
all.muts[, subject := sample.to.subject.map[sample]]
all.muts[, id := paste(chr, pos, refnt, altnt)]

final.muts <- c()

# Filtering for things like recurrence/clustering is done separately
# for SNVs and indels.  Perhaps clustering would be better if the two
# were combined?
for (mt in muttypes) {
    # Determine confidence class for every SNV
    muts <- all.muts[muttype == mt]

    cat('Initial statistics:', sum(muts$pass),
        'passed', sum(muts$rescue), 'rescued\n')
    

    # Remove mutations called in 2 different subjects.
    # This does NOT filter SNVs called more than once in the same
    # individual (=likely lineage marker), in which case one of the
    # multply-called mutations is retained. That filter is applied later.
    cat('Raw recurrence rates:\n')
    print(table(table(muts$id)))
    
    cat('Recurrence x subject table\n')
    z <- split(muts$subject, muts$id)
    subjects <- sapply(z, function(v) length(unique(v)))
    recs <- sapply(z, length)
    print(addmargins(table(recs, subjects)))
    muts[, rec.filter := subjects[id] > 1]


    # nearby points created by a single sample are more likely to
    # be artifacts. Remove the whole cluster, because it is often
    # true that the entire cluster is caused by the same few reads
    # that probably align poorly or are clipped.
    filter.single.sample.clusters <- function(muts, threshold=300) {
        by.sample <- split(muts, muts$sample)
        muts <- do.call(rbind, lapply(by.sample, function(d) {
            do.call(rbind, lapply(split(d, d$chr), function(dchr) {
                dchr <- dchr[order(dchr$pos),]
                up <- c(Inf, diff(dchr$pos))
                down <- c(diff(dchr$pos), Inf)
                dchr$nearest <- pmin(up, down)
                dchr
            }))
        }))
        # return back to original order
        muts <- muts[order(chr, pos)]
        muts[, cluster.filter := nearest < threshold]
        cat(sprintf('Removing %d sites within %d bp in the same sample\n',
            sum(muts$cluster.filter), threshold))
        print(addmargins(table(muts$sample, muts$cluster.filter)))
        muts
    }
    
    muts <- filter.single.sample.clusters(muts, threshold=args$cluster_filter_bp)

    # if there are any recurrent mutations remaining,
    # then they recur in the same subject and are likely
    # lineage-related true mutations. it's important not to count them as
    # multiple independent mutations or high clonality lineage markers
    # could drive enrichment signals.
    muts[, lineage.filter := duplicated(id)]
    muts[, final.filter := rec.filter | cluster.filter | lineage.filter]

    final.muts <- rbind(final.muts, muts[pass == TRUE & rescue == TRUE & final.filter == TRUE])

    if (!is.null(args$separate_files)) {
        for (passtype in c('pass', 'pass_and_rescue')) {
            outfile=output.tables[paste0(mt, '_', passtype)]
            if (file.exists(outfile))
                stop(paste('output file', outfile, 'already exists, please delete it first'))
            if (passtype == 'pass')
                outmuts <- muts[pass == TRUE]
            if (passtype == 'pass_and_rescue')
                outmuts <- muts[pass == TRUE | rescue == TRUE]

            cat('writing', outfile, '\n')
            data.table::fwrite(outmuts[final.filter == FALSE], file=outfile)
        }
    }
}

data.table::fwrite(all.muts, file=main.output.table)
