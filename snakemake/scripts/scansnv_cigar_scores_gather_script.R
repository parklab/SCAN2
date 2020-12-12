library(scan2)

gt <- do.call(concat, lapply(snakemake@input[['rdas']], function(f) {
    print(f)
    get(load(f))
}))

save(gt, file=snakemake@output[[1]])
