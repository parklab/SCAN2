#!/usr/bin/env Rscript

# same as torda.R, but only analyzes the X chromosome

require('scansnv')
args <- commandArgs(TRUE)

if (length(args) != 2) {
    cat("outputformatprefix is a printf-style string that must contain\n")
    cat("a %d placeholder.  the %d will be replaced by a chromosome\n")
    cat("integer for the chromosome-specific files and 'all' for the\n")
    cat("full RData file.  since both .bin (per-chrom) and .rda (full)\n")
    cat("formats are written, the prefix should not include a file\n")
    cat("extension\n")
    stop("Rscript torda.R training_data.txt outputformatprefix")
}

f <- args[1]
outfmt <- args[2]

if (length(grep("%s", outfmt)) == 0)
    stop("outputformat must contain an instance of %d")

cat(sprintf("reading data from %s...\n", f))
data <- read.table(f, header=TRUE, stringsAsFactors=FALSE)
cat(sprintf("%d rows, summaries:\n", nrow(data)))
diffs <- diff(data$pos)
do.call(rbind, list(
    hap1=summary(data$hap1),
    hap2=summary(data$hap2),
    dp=summary(data$dp),
    `distance between adjacent sites`= summary(diffs)))
cat("table of left vs. right SNP alleles")
table(data$phgt)

outf <- sprintf("%s.rda", sub("%s", "all", outfmt))
if (file.exists(outf))
    stop(sprintf("output file %s already exists, aborting", outf))
cat(sprintf("writing full data (%d rows) to %s\n", nrow(data), outf))
save(data, file=outf)

for (chr in "X") {
    outf <- sprintf("%s.bin", sprintf(outfmt, chr))
    if (file.exists(outf))
        stop(sprintf("output file %s already exists, aborting", outf))

    dchr <- data[data$chr==chr,]
    cat(sprintf("saving %d rows to %s...\n", nrow(dchr), outf))
    write.bin(dchr, fname=outf)
    #save(dchr, file=outf)
}
