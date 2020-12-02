library(viridis)

# Sum over all chunked tables
dptab <- Reduce('+', lapply(snakemake@input[['rdas']], function(infile) {
    load(infile)
    dptab
}))


# (0, 0) is usually very large, but largely reflects Ns and
# unalignable regions.
dptab.tmp <- dptab
dptab.tmp[1,1] <- 0

pdf(snakemake@output[['pdf']])
filled.contour(x=0:200, y=0:200, dptab.tmp,
    color.palette=viridis,
    xlab=snakemake@wildcards[['sample']],
    ylab=snakemake@config[['bulk_sample']],
    main='Sequencing depth')
dev.off()

save(dptab, file=snakemake@output[['rda']])
