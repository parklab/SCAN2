# save data to a variable named by this string
varname <- snakemake@params[['varname']]

#set to above
tmp <- do.call(rbind,
    lapply(snakemake@input, function(f) {
        load(f)
        gt$somatic
    })
)
# equivalent to `varname` <- tmp, where `varname` substitutes the string value
assign(varname, tmp)

save(list=varname, file=snakemake@output[[1]])
