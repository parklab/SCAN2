library(scan2) # provides resample.hsnps()

vcf <- read.table(snakemake@input[['vcf']],
    stringsAsFactors=TRUE, header=F,
    comment="#", colClasses=c(V1='character'))
colnames(vcf)[c(1:5)] <- c('chr', 'pos', 'dbsnp', 'refnt', 'altnt')
vcf <- vcf[vcf$chr == snakemake@wildcards[['chr']], c(1:2,4:5)]

som <- read.table(snakemake@input[['somatic_pos']], header=T,
    stringsAsFactors=FALSE, colClasses=c(chr='character'))

str(vcf)
str(som)
print(snakemake@wildcards[['chr']])
print(snakemake@config[['resample_M']])

rs <- resample.hsnps(som, vcf,
    chrom=snakemake@wildcards[['chr']], M=snakemake@config[['resample_M']])

pdf(snakemake@output[['pdf']])
plot(rs$dist.s$mids, rs$dist.s$density,
    type='l', col='red', lwd=2, xlab='log10(dist between sites)',
    ylab='Density', main='Intra-hSNP distances; hSNP-sSNV distances')
lines(rs$dist.h$mids, rs$dist.h$density, col=1, lwd=2)
dev.off()
hsnp.sample <- vcf[rs$selection$keep,]
save(hsnp.sample, file=snakemake@output[['rda']])
write.table(hsnp.sample,
    file=snakemake@output[['tab']], quote=F, row.names=FALSE, sep='\t')
save(rs, file=snakemake@output[['resample_rda']])
