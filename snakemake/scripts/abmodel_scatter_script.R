library(scan2)
alim <- c(-7, 2)
blim <- c(2, 4)
clim <- c(-7, 2)
dlim <- c(2, 6)
if (snakemake@params[['step']] > 1) {
    load(snakemake@params[['paramfile']])
    alim <- bounds[,1]
    blim <- bounds[,2]
    clim <- bounds[,3]
    dlim <- bounds[,4]
}
load(snakemake@input[['training']])
data <- data[data$chr == snakemake@params[['chr']],]
ctx <- abmodel.approx.ctx(x=data$pos,
    y=data$hap1, d=data$hap1+data$hap2,
    hsnp.chunksize=snakemake@config[['abmodel_hsnp_tile_size']]
)

logp.samples <- abmodel.sample(
    n=snakemake@config[['abmodel_samples_per_chunk']],
    alim=alim, blim=blim, clim=clim, dlim=dlim,
    ctx=ctx,
    seed=snakemake@params[['seed']])
save(logp.samples, file=snakemake@output[[1]])
