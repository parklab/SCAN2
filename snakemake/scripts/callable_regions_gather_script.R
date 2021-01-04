library(viridis)

# Get the clamp DP (all DP values > clamp are set to clamp)
load(snakemake@input[['rdas']][1]) # loads clamp.dp, dptab

# Sum over all chunked tables (overwrites dptab from above)
dptab <- Reduce('+', lapply(snakemake@input[['rdas']], function(infile) {
    load(infile)
    dptab
}))


# (0, 0) is usually huge, but largely reflects Ns and
# unalignable regions.
dptab.tmp <- dptab
dptab.tmp[1,1] <- 0

pdf(snakemake@output[['pdf']])
filled.contour(x=0:100, y=0:100, dptab.tmp[1:101,1:101],
    color.palette=viridis,
    xlab=snakemake@wildcards[['sample']],
    ylab=snakemake@config[['bulk_sample']],
    main='Sequencing depth')
dev.off()

save(dptab, file=snakemake@output[['rda']])
