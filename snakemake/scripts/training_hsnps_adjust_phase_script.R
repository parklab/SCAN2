library(scan2)
library(data.table)

unadjusted <- fread(snakemake@input[['tab']])

if (file.exists(snakemake@output[['tab']]))
    stop(paste("output file", snakemake@output[['tab']],
        "already exists, please delete it first"))

if (snakemake@config[['parsimony_phasing']]) {
    unadjusted$af <- unadjusted$hap1 / unadjusted$dp
    adjusted <- adjust.phase(unadjusted)
} else {
    adjusted <- unadjusted
}

fwrite(adjusted, file=snakemake@output[['tab']], sep='\t')
