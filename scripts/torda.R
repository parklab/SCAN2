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

if (length(grep("%d", outfmt)) == 0)
    stop("outputformat must contain an instance of %d")

# No longer in separate file. Only used here.

# write out the position, hap1, dp info for consumption by the C
# mkl-gridfit program.  format:
# N <integer> number of rows
# pos N*<integer> vector of positions
# hap1 N*<integer> vector of hap1s
# dp N*<integer> vector of DPs
write.bin <- function(df, fname) {
    f <- file(fname, 'wb')
    writeBin(nrow(df), f)
    writeBin(df$pos, f)
    writeBin(df$hap1, f)
    writeBin(df$hap1 + df$hap2, f)
    close(f)
}

cat(sprintf("reading data from %s...\n", f))
data <- read.table(f, header=TRUE, stringsAsFactors=FALSE)
cat(sprintf("%d rows, summaries:\n", nrow(data)))
diffs <- unlist(lapply(1:22, function(chr) diff(data[data$chr==chr,]$pos)))
do.call(rbind, list(
    hap1=summary(data$hap1),
    hap2=summary(data$hap2),
    dp=summary(data$dp),
    `distance between adjacent sites`= summary(diffs)))
cat("table of left vs. right SNP alleles")
table(data$phgt)

outf <- sprintf("%s.rda", sub("%d", "all", outfmt))
if (file.exists(outf))
    stop(sprintf("output file %s already exists, aborting", outf))
cat(sprintf("writing full data (%d rows) to %s\n", nrow(data), outf))
save(data, file=outf)

for (chr in 1:22) {
    outf <- sprintf("%s.bin", sprintf(outfmt, chr))
    if (file.exists(outf))
        stop(sprintf("output file %s already exists, aborting", outf))

    dchr <- data[data$chr==chr,]
    cat(sprintf("saving %d rows to %s...\n", nrow(dchr), outf))
    write.bin(dchr, fname=outf)
    #save(dchr, file=outf)
}
