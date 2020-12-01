x <- do.call(rbind, 
    lapply(snakemake@input, function(f) {
        load(f)
        logp.samples
    }))
dn <- dimnames(x)
x <- as.matrix(x)
# swapping columns (1,2) and (3,4) to force b < d.
# ifelse returns a value the same shape as the first argument
logi.mat <- matrix(rep(x[,2] < x[,4], times=5), ncol=5)
x <- as.data.frame(ifelse(logi.mat, x, x[,c(3,4,1,2,5)]))
dimnames(x) <- dn
x <- x[order(x[,5], decreasing=TRUE),]

# The highest logp value (x[,5]) is the best fit
fit <- x[1,,drop=FALSE]
chr <- snakemake@wildcards[['chr']]
save(chr, fit, file=snakemake@output[['fit']])

# Use the top 50 logp values to build a new parameter range
x[,2] <- log10(x[,2])   # b and d bounds are in log10 space
x[,4] <- log10(x[,4])
bounds <- apply(head(x[,-5], 50), 2, range)
colnames(bounds) <- colnames(x)[-5]
save(bounds, file=snakemake@output[['range']])
