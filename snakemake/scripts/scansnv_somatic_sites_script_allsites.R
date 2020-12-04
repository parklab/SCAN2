# Use a very sensitive definition of somatic site here
# At least: 0 alt reads in bulk, not a no-call, not in dbsnp
            # and at least 2 non-bulk reads.
tab <- read.table(snakemake@input[[1]], stringsAsFactors=FALSE, header=TRUE,
                comment="#", colClasses=c(chr='character'))
bulk.sample <- make.names(snakemake@config[['bulk_sample']])
bulk.idx <- which(colnames(tab) == bulk.sample)
bulk.alt <- bulk.idx +  2
sc.alts <- which(grepl("alt", colnames(tab)) &
                 colnames(tab) != colnames(tab)[bulk.alt] &
                 colnames(tab) != "altnt")
cat("Using data:\n")
for (i in 1:ncol(tab)) {
    s <- ifelse(i %in% sc.alts,
        "[SC]",
        ifelse(i == bulk.alt, "[BULK]", ""))
    cat(sprintf("%8s %s\n", s, colnames(tab)[i]))
}

# in allsites mode, calculate cigars/AB/etc for all positions
candidate.somatics <- tab

write.table(candidate.somatics[,c(1:2, 4:5)], sep="\t",
    quote=FALSE, row.names=FALSE, file=snakemake@output[[1]])
