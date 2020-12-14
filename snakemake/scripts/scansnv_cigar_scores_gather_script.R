library(scan2)

gt <- do.call(concat, lapply(snakemake@input[['rdas']], function(f) {
    print(f)
    get(load(f))
}))

gt <- add.static.filters(gt)
gt <- compute.fdr.priors(gt)
gt <- compute.fdr(gt)

save(gt, file=snakemake@output[[1]])
