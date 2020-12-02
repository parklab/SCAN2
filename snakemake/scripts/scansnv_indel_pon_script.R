load(snakemake@config[['somatic_indel_pon']])   # loads 'pon'
load(snakemake@input[['rda']])                  # loads 'somatic'

# PON has multiallelic sites, so cannot include refnt/altnt in merge criteria.
# XXX: maybe refnt, though?
somatic <- merge(somatic, pon[,-(3:5)], by=c('chr','pos'))
somatic$pass.prepon <- somatic$pass
somatic$pass <- somatic$pass & somatic$dp >= 10 &
    (somatic$unique.donors <= 1 | somatic$max.out < 2)
save(somatic, file=snakemake@output[['rda']])
