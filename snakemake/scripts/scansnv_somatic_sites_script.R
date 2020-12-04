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

if (snakemake@wildcards[['muttype']] == "snv") {
    candidate.somatics <-
        tab[tab[,bulk.alt] <= snakemake@config[['max_bulk_alt']] &
            tab[,bulk.idx] == '0/0' &
            tab$dbsnp == '.' &
            rowSums(as.matrix(tab[,sc.alts])) >= snakemake@config[['min_sc_alt']] &
            tab$chr == snakemake@wildcards[['chr']],]
} else if (snakemake@wildcards[['muttype']] == "mosaic_snv") {
    bulk.vaf <- tab[,bulk.alt] / (tab[,bulk.alt] + tab[,bulk.alt-1])
    candidate.somatics <-
        tab[tab[,bulk.alt] > 0 &
            bulk.vaf < snakemake@config[['max_bulk_mosaic_vaf']] &
            tab$dbsnp == '.' &
            rowSums(as.matrix(tab[,sc.alts])) >= snakemake@config[['min_sc_alt']] &
            # To help remove germline: at least one cell must have 0 alt reads
            rowSums(as.matrix(tab[,sc.alts]) == 0) > 0 &
            tab$chr == snakemake@wildcards[['chr']],]
} else
    stop(sprintf("params.muttype=%s is unrecognized. must be 'snv' or 'mosaic_snv'", snakemake@wildcards[['muttype']]))

write.table(candidate.somatics[,c(1:2, 4:5)], sep="\t",
    quote=FALSE, row.names=FALSE, file=snakemake@output[[1]])
