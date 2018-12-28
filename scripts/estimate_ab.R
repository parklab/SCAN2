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
    # We're predicting AB at hSNPs as a control. These hSNPs are derived
    # from the training data, so we need to exclude them in this computation.
    # N.B. this isn't a perfect control scenario: the hSNPs were still used
    # to fit the correlation function parameters; however, as long as a
    # small fraction of hSNPs are used here, the fit would not have been
    # affected noticeably anyway.
    if (analysis.type == 'hsnp_spikein') {
str(hsnps)
str(sites)
        exclusion <- paste(hsnps$chr, hsnps$pos, hsnps$altnt) %in%
                     paste(sites$chr, sites$pos, sites$altnt)
        hsnps <- hsnps[!exclusion,]
        cat(sprintf("hsnp_spikein analysis: withheld %d sites on chromosome %s\n",
            sum(exclusion), chrom))
    }
    fit <- fits[[chrom]]

    cat(sprintf("inferring AB for %d sites on chr%s:%d-%d\n", 
        nrow(sites), chrom, min(sites$pos), max(sites$pos)))
    system.time(z <- infer.gp(ssnvs=sites, fit=fit,
        hsnps=hsnps, chunk=1, flank=1e5))
    cbind(sites, z)
}))

save(ab, file=outfile)
