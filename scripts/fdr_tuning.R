#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 11)
    stop("usage: mmq60.tab mmq1.tab training.rda bulk_sample_name sc_sample_name output.rda { somatic | spikein} min_sc_alt min_sc_dp min_bulk_dp somatic_sites1 [ somatic_sites2 ... somatic_sitesN ]\n    note: in spikein mode, sites in somatic_sites1 are not filtered due to alt reads detected in bulk.")

library(scan2)

mmq60.file <- args[1]
mmq1.file <- args[2]
training.file <- args[3]
bulk.sample <- make.names(args[4])
sc.sample <- make.names(args[5])
outfile <- args[6]
spikein <- args[7]
min.sc.alt <- as.integer(args[8])
min.sc.dp <- as.integer(args[9])
min.bulk.dp <- as.integer(args[10])
site.files <- args[-(1:10)]

if (spikein != 'somatic' & spikein != 'spikein')
    stop(sprintf("expected 'somatic' or 'spikein' but got %s", spikein))

spikein_mode = spikein == 'spikein'

somatic.sites <- do.call(rbind, lapply(site.files,
    function(f)
        read.table(f, header=T, stringsAsFactors=FALSE,
            colClasses=c(chr='character'))
))


nrow(somatic.sites)
if (spikein_mode) {
    protected <- read.table(site.files[1], header=T, stringsAsFactors=F, colClasses=c(chr='character'))
nrow(protected)
}


load(training.file)  # loads 'data'

# Step 1: read a few rows just to get column names
hmq <- read.table(mmq60.file, header=T, stringsAsFactors=F,
    colClasses=c(chr='character'), nrow=10)
tot.cols <- ncol(hmq)
sc.idx <- which(colnames(hmq) == sc.sample)
bulk.idx <- which(colnames(hmq) == bulk.sample)
cat("Using data:\n")
for (i in 1:ncol(hmq)) {
    s <- ifelse(i == sc.idx, "[SC]", ifelse(i == bulk.idx, "[BULK]", ""))
    cat(sprintf("%8s %s\n", s, colnames(hmq)[i]))
}

# Step 2: now read only the NEEDED columns. This is critical for
# very large data sets where these tables can easily require 20G of
# RAM or more to store in R.
# N.B. "NULL" means read.table will skip the column, NA means it
# will try to automatically determine the column's data type.
cols.to.read <- rep("NULL", tot.cols)
# First 7 are chr, pos, dbsnp, refnt, altnt, mq, mqrs
cols.to.read[1:7] <- c('character', 'integer', rep('character', 3), 'numeric', 'numeric')
# Read 3 columns for the single cell, 3 columns for bulk
cols.to.read[sc.idx + 0:2] <- c('character', 'integer', 'integer')
cols.to.read[bulk.idx + 0:2] <- c('character', 'integer', 'integer')
cat(sprintf("reading %d columns from %s..\n",
    sum(is.na(cols.to.read) | cols.to.read != 'NULL'), mmq60.file))
hmq <- read.table(mmq60.file, header=T, stringsAsFactors=F,
    colClasses=cols.to.read)
new.sc.idx <- which(colnames(hmq) == sc.sample)
new.bulk.idx <- which(colnames(hmq) == bulk.sample)
str(hmq)

hmq$dp <- hmq[,new.sc.idx+1] + hmq[,new.sc.idx+2]
hmq$af <- hmq[,new.sc.idx+2] / hmq$dp
hmq$bulk.dp <- hmq[,new.bulk.idx+1] + hmq[,new.bulk.idx+2]


# The mmq=1 table is even larger (sometimes 2-3x larger) than
# the mmq=1 table. The only data we need from mmq=1 is the number
# of alt reads in bulk.
cols.to.read <- rep("NULL", tot.cols)
cols.to.read[1:7] <- c('character', 'integer', rep('character', 3), 'numeric', 'numeric')
# Read 3 columns for the single cell, 3 columns for bulk
cols.to.read[bulk.idx + 2] <- 'integer'
cat(sprintf("reading %d columns from %s..\n",
    sum(is.na(cols.to.read) | cols.to.read != 'NULL'), mmq1.file))
lmq <- read.table(mmq1.file, header=T, stringsAsFactors=F,
    colClasses=cols.to.read)
colnames(lmq)[8] <- 'lowmq.bulk.alt'
str(lmq)


# In joint calling mode, there can be an unbounded number of somatic
# candidates from other cells.  The majority of candidates are WGA
# errors that occurred in other cells, which will have VAF=0 in this
# cell.
# It is important for the FDR tuning procedure to only consider
# candidates for *this* sample.
# XXX: unfortunate sloppy coding: all of this must match the candidate
# selection logic in r-scan2::genotype_somatic().
somatic.candidates <- merge(somatic.sites, hmq, all.x=T)
somatic.candidates <- merge(somatic.candidates, lmq[,c('chr', 'pos', 'lowmq.bulk.alt')], all.x=T)

# do not apply the lowmq alt test to spikeins, since these should have
# supporting reads in bulk.
somatic.candidates$spikein.exception <- FALSE
if (spikein_mode)
    somatic.candidates$spikein.exception <-
        paste(somatic.candidates$chr, somatic.candidates$pos) %in%
        paste(protected$chr, protected$pos)
table(somatic.candidates$spikein.exception)
str(somatic.candidates)
table(somatic.candidates$spikein.exception |
                        is.na(somatic.candidates$lowmq.bulk.alt) |
                        somatic.candidates$lowmq.bulk.alt == 0)

somatic.candidates <-
    somatic.candidates[somatic.candidates[,new.sc.idx+2] >= min.sc.alt &
                       somatic.candidates$dp >= min.sc.dp &
                       somatic.candidates$bulk.dp >= min.bulk.dp &
                       (somatic.candidates$spikein.exception |
                        is.na(somatic.candidates$lowmq.bulk.alt) |
                        somatic.candidates$lowmq.bulk.alt == 0),]

hsnps <- merge(data[,c('chr', 'pos')], hmq, all.x=T)
hsnps <- hsnps[hsnps[,new.sc.idx+2] >= min.sc.alt,]

print("SOMATIC  CANDIDATES")
str(somatic.candidates)
print("HSNPS")
str(hsnps)
fdr.tuning <- get.fdr.tuning.parameters(somatic.candidates, hsnps)
save(fdr.tuning, file=outfile)
