#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5)
    stop("usage: mmq60.tab training.rda sc_sample_name output.rda somatic_sites1 [ somatic_sites2 ... somatic_sitesN ]")

library(scansnv)

mmq60.file <- args[1]
training.file <- args[2]
sc.sample <- make.names(args[3])
outfile <- args[4]
site.files <- args[-(1:4)]


somatic.sites <- do.call(rbind, lapply(site.files,
    function(f)
        read.table(f, header=T, stringsAsFactors=FALSE,
            colClasses=c(chr='character'))
))

load(training.file)  # loads 'data'

hmq <- read.table(mmq60.file, header=T, stringsAsFactors=T,
    colClasses=c(chr='character'))

sc.idx <- which(colnames(hmq) == sc.sample)
cat("Using data:\n")
for (i in 1:ncol(hmq)) {
    s <- ifelse(i == sc.idx, "[SC]", "")
    cat(sprintf("%8s %s\n", s, colnames(hmq)[i]))
}
hmq$dp <- hmq[,sc.idx+1] + hmq[,sc.idx+2]
hmq$af <- hmq[,sc.idx+2] / hmq$dp

somatic.candidates <- merge(somatic.sites, hmq, all.x=T)
hsnps <- merge(data[,c('chr', 'pos')], hmq, all.x=T)

fdr.tuning <- get.fdr.tuning.parameters(somatic.candidates, hsnps)
save(fdr.tuning, file=outfile)
