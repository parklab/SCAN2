#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 11)
    stop("usage: mmq60.tab mmq1.tab training.rda bulk_sample_name sc_sample_name output.rda { somatic | spikein} min_sc_alt min_sc_dp min_bulk_dp somatic_sites1 [ somatic_sites2 ... somatic_sitesN ]\n    note: in spikein mode, sites in somatic_sites1 are not filtered due to alt reads detected in bulk.")

library(scansnv)

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
if (spikein_mode)
{
    protected <- read.table(site.files[1], header=T, stringsAsFactors=F, colClasses=c(chr='character'))
nrow(protected)
}


load(training.file)  # loads 'data'

hmq <- read.table(mmq60.file, header=T, stringsAsFactors=F,
    colClasses=c(chr='character'))

sc.idx <- which(colnames(hmq) == sc.sample)
bulk.idx <- which(colnames(hmq) == bulk.sample)
cat("Using data:\n")
for (i in 1:ncol(hmq)) {
    s <- ifelse(i == sc.idx, "[SC]", ifelse(i == bulk.idx, "[BULK]", ""))
    cat(sprintf("%8s %s\n", s, colnames(hmq)[i]))
}
hmq$dp <- hmq[,sc.idx+1] + hmq[,sc.idx+2]
hmq$af <- hmq[,sc.idx+2] / hmq$dp
hmq$bulk.dp <- hmq[,bulk.idx+1] + hmq[,bulk.idx+2]

lmq <- read.table(mmq1.file, header=T, stringsAsFactors=F,
    colClasses=c(chr='character'))

# In joint calling mode, there can be an unbounded number of somatic
# candidates from other cells.  The majority of candidates are WGA
# errors that occurred in other cells, which will have VAF=0 in this
# cell.
# It is important for the FDR tuning procedure to only consider
# candidates for *this* sample.
# XXX: unfortunate sloppy coding: all of this must match the candidate
# selection logic in r-scansnv::genotype_somatic().
somatic.candidates <- merge(somatic.sites, hmq, all.x=T)
colnames(lmq)[bulk.idx+2] <- 'lowmq.bulk.alt'
nrow(somatic.candidates)
somatic.candidates <- merge(somatic.candidates, lmq[,c('chr', 'pos', 'lowmq.bulk.alt')], all.x=T)

# do not apply the lowmq alt test to spikeins, since these should have
# supporting reads in bulk.
nrow(somatic.candidates)
somatic.candidates$spikein.exception <- FALSE
nrow(somatic.candidates)
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
    somatic.candidates[somatic.candidates[,sc.idx+2] >= min.sc.alt &
                       somatic.candidates$dp >= min.sc.dp &
                       somatic.candidates$bulk.dp >= min.bulk.dp &
                       (somatic.candidates$spikein.exception |
                        is.na(somatic.candidates$lowmq.bulk.alt) |
                        somatic.candidates$lowmq.bulk.alt == 0),]

hsnps <- merge(data[,c('chr', 'pos')], hmq, all.x=T)
hsnps <- hsnps[hsnps[,sc.idx+2] >= min.sc.alt,]

print("SOMATIC  CANDIDATES")
str(somatic.candidates)
print("HSNPS")
str(hsnps)
fdr.tuning <- get.fdr.tuning.parameters(somatic.candidates, hsnps)
save(fdr.tuning, file=outfile)
