library(scan2)

unadjusted <- read.table(snakemake@input[['tab']],
    header=TRUE, stringsAsFactors=FALSE)

if (file.exists(snakemake@output[['tab']]))
    stop(paste("output file", snakemake@output[['tab']],
        "already exists, please delete it first"))

unadjusted$af <- unadjusted$hap1 / unadjusted$dp
adjusted <- adjust.phase(unadjusted)

write.table(adjusted, file=snakemake@output[['tab']],
    quote=F, row.names=FALSE, sep='\t')
