#!/usr/bin/env Rscript

require('scansnv')
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5)
    stop("usage: somatic_ab.R fits.rda training.rda positions.txt { hsnp_spikein | somatic } output.rda\n")

fitfile <- args[1]
trainfile <- args[2]
posfile <- args[3]
analysis.type <- args[4]
outfile <- args[5]

if (file.exists(outfile))
    stop(sprintf("output file %s already exists, please delete it first", outfile))

if (analysis.type != "hsnp_spikein" & analysis.type != "somatic")
    stop("analysis type (argument 4) must be either 'hsnp_spikein' or 'somatic'")

load(fitfile)    # loads 'fits'
load(trainfile)  # loads 'data'
sites <- read.table(posfile, sep="\t", header=TRUE,
    stringsAsFactors=FALSE, colClasses=c(chr='character'))



# This is not for parallelization. Each chromosome is fit with a different
# correlation function, so the correct parameters must be supplied.
ab <- do.call(rbind, lapply(unique(sites$chr), function(chrom) {
    hsnps <- data[data$chr == chrom,]
    fit <- fits[[chrom]]

    cat(sprintf("inferring AB for %d sites on chr%s:%d-%d\n", 
        nrow(sites), chrom, min(sites$pos), max(sites$pos)))
    system.time(z <- infer.gp(ssnvs=sites, fit=fit,
        hsnps=hsnps, chunk=1, flank=1e5, verbose=FALSE,
        spikein=analysis.type == 'hsnp_spikein'))
    cbind(sites, z)
}))

save(ab, file=outfile)
