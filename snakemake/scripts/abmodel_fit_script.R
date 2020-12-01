x <- lapply(snakemake@input, function(f) {
    load(f)
    list(chr=chr, fit=fit)
})
fits <- lapply(x, function(xx) xx$fit)
names(fits) <- lapply(x, function(xx) xx$chr)
save(fits, file=snakemake@output[[1]])
