mosaic <- do.call(rbind,
    lapply(snakemake@input[['filelist']], function(f) {
        load(f)
        gt$somatic
    })
)

save(mosaic, file=snakemake@output[[1]])
