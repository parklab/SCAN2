#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser()

parser$add_argument('output_prefix', type='character', metavar='PREFIX',
    help='This script creates 4 tables: one for each combination of pass vs. rescue and snv vs. indel. These tables will be named PREFIX_{pass|pass_and_rescue}_{snv|indel}.txt')
parser$add_argument('--metadata', metavar='FILE', default=NULL,
    help='CSV file mapping single cell sample IDs (as in the sample column of --muts files) to individual IDs (i.e., brain donors). If this file is not specified, then each cell will be treated as though it comes from a different donor. This affects the recurrence filter, which removes any mutation occurring in more than one individual. If a mutation occurs in multiple cells from the same individual, such as a lineage marker, then one of them is retained but this can only be determined if sample->subject metadata is given.')
parser$add_argument('--muts', action='append', metavar='FILE',
    help='CSV file of somatic mutations with at least the following columns: sample, chr, pos, refnt, altnt, muttype, pass, rescue. This argument can be specified multiple times to combine tables from multiple runs. Must be specified at least once.')
parser$add_argument("--cluster-filter-bp", metavar='INT', type='integer', default=50,
    help='Remove mutations that are within INT base pairs of the nearest mutation of the same type (i.e., snv or indel) AND in the same sample.')

args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

if (length(args$muts) < 1)
    stop('at least one --muts argument must be given')


library(scan2)

step4.exact_recurrence_filter=F      # filter out SNVs called in more than one donor.
step5.clustered_mut_filter=F         # remove SNVs within 50bp in the same sample
                                     # these are mostly dinuc subs.
step6.finalize=F    
step7.save_output=T                  # write (but don't overwrite) the rda file

# XXX: NEED METADATA FOR SAMPLE -> DONOR MAPPING

# Filtering for things like recurrence/clustering is done separately
# for SNVs and indels.  Perhaps clustering would be better if the two
# were combined?

all.muts <- do.call(rbind, lapply(args$muts, function(mutfile)
    data.table::fread(mutfile)[, .(sample, chr, pos, refnt, altnt, muttype, pass, rescue)]))

# Default (no metadata): each cell is treated as coming from a unique subject
sample.to.subject.map <- setNames(unique(all.muts$sample), unique(all.muts$sample))
if (!is.null(args$metadata)) {
    meta <- fread(args$metadata)
    sample.to.subject.map <- setNames(meta$sample, meta$subject)
}


all.muts <- all.muts[order(chr, pos),]
all.muts[, subject := sample.to.subject.map[sample]]
all.muts[, id := paste(chr, pos, refnt, altnt)]

muttypes <- c('snv', 'indel')
for (mt in muttypes) {
    # Determine confidence class for every SNV
    muts <- all.muts[muttype == mt]

    cat('Initial statistics:', sum(muts$pass),
        'passed', sum(muts$rescue), 'rescued\n')
    

    # Remove mutations called in 2 different subjects.
    # This does NOT filter SNVs called more than once in the same
    # individual (=likely lineage marker), in which case one of the
    # multply-called mutations is retained. That filter is applied later.
    cat('Raw recurrence rates (PTA):\n')
    print(table(table(nmut$id[nmut$passA | nmut$passB])))
    
    cat('Recurrence x donor table (all muts)\n')
    z <- split(nmut$donor, nmut$id)
    donors <- sapply(z, function(v) length(unique(v)))
    recs <- sapply(z, length)
    print(addmargins(table(recs,donors)))
    nmut$rec.filter <- donors[nmut$id] > 1

    cat('Recurrence x donor table (any passA,B,M)\n')
    z <- split(nmut$donor[nmut$passA | nmut$passB],
            nmut$id[nmut$passA | nmut$passB])
    donors <- sapply(z, function(v) length(unique(v)))
    recs <- sapply(z, length)
    print(addmargins(table(recs,donors)))


# nearby points created by a single sample are more likely to
# be artifacts. Remove the whole cluster, because it is often
# true that the entire cluster is caused by the same few reads
# that probably align poorly or are clipped.
if (step5.clustered_mut_filter) {
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
        muts <- muts[order(muts$chr, muts$pos),]
        filt.name <- paste0('clustered.filt.', threshold)
        muts[[filt.name]] <- muts$nearest < threshold
        cat(sprintf('Removing %d sites within %d bp in the same sample\n',
            sum(muts[[filt.name]]), threshold))
        print(addmargins(table(muts$sample, muts[[filt.name]])))
        muts
    }
    
    nmut <- filter.single.sample.clusters(nmut, threshold=50)
}


if (step6.finalize) {
    nmut <- nmut[nmut$passA | nmut$passB | nmut$passM,]

    # if there are any recurrences left, be sure to only pick one. if the
    # mutation is observed in both MDA and PTA, prefer the PTA one. otherwise,
    # keep one of any copy.
    # these recurrent mutations that occur in the same subject are likely
    # lineage-related true mutations. it's important not to count them as
    # multiple independent mutations or widely-inherited lineage markers, e.g.,
    # could drive enrichment signals.
    ord.samples <- c(names(ss.pta), names(ss.mda))
    nmut <- do.call(rbind, lapply(ord.samples, function(sn)
        nmut[nmut$sample == sn,]))
    nmut$lineage.filter <- duplicated(nmut$id)
    #nmut <- nmut[!nmut$duplicated,]
    nmut <- nmut[order(nmut$chr, nmut$pos),]
    nmut$final.filter <- nmut$rec.filter | nmut$clustered.filt.50 | nmut$lineage.filter

    # give a pass status that assigns each mutation to its most confident set
    # the priority of pass status is A > B > M
    # so if a mutation is both passA and passB, it gets assigned to the passA set.
    nmut$status <- ifelse(nmut$passA, 'A',
        ifelse(nmut$passB, 'B',
                    ifelse(nmut$passM, 'M', 'ERROR')))
}

if (step7.save_output) {
    cat('step7\n')
    rdafile='/n/data1/hms/dbmi/park/jluquette/pta/integrated_calls/calls_20220506_clus50bp_rec1donor_passABM_all_sites_in_table.rda'
    if (file.exists(rdafile))
        stop(paste('output RDA file', rdafile, 'already exists, please delete it first'))
    save(nmut, file=rdafile)

    csvfile='/n/data1/hms/dbmi/park/jluquette/pta/integrated_calls/calls_20220506_clus50bp_rec1donor_passABM_all_sites_in_table.csv'
    if (file.exists(csvfile))
        stop(paste('output csv file', csvfile, 'already exists, please delete it first'))
    fwrite(nmut, file=csvfile)
}
