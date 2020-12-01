cigars <- do.call(rbind, lapply(snakemake@input, function(f) 
    read.table(f, header=T)))
write.table(cigars, file=snakemake@output[[1]],
    quote=FALSE, row.names=FALSE, sep="\t")
