library(scan2)

gt <- do.call(concat, lapply(snakemake@input[['rdas']], function(f) {
    print(f)
    get(load(f))
}))

gt <- add.ab.fits(gt, snakemake@input[['ab_fits']])
gt <- add.training.data(gt, snakemake@input[['ab_training']])
gt <- resample.training.data(gt, M=snakemake@config[['resample_M']])

save(gt, file=snakemake@output[[1]])
