tab <- read.table(snakemake@input[[1]],
    stringsAsFactors=FALSE, header=TRUE,
    comment="#", colClasses=c(chr='character'))
bulk.sample <- make.names(snakemake@config[['bulk_sample']])
bulk.idx <- which(colnames(tab) == bulk.sample)
bulk.alt <- bulk.idx +  2
sc.alts <- which(grepl("alt", colnames(tab)) &
                 colnames(tab) != colnames(tab)[bulk.alt] &
                 colnames(tab) != "altnt")
cat("Using data:\n")
for (i in 1:ncol(tab)) {
    s <- ifelse(i %in% sc.alts, "[SC]", ifelse(i == bulk.alt, "[BULK]", ""))
    cat(sprintf("%8s %s\n", s, colnames(tab)[i]))
}

bulk.vaf <- tab[,bulk.alt] / (tab[,bulk.alt] + tab[,bulk.alt-1])
candidate.mosaics <-
    tab[tab[,bulk.alt] > 0 &
        bulk.vaf < 0.3 &
        rowSums(as.matrix(tab[,sc.alts])) >= snakemake@config[['min_sc_alt']] &
        # To help remove germline: at least one cell must have 0 alt reads,
        # which operates as a quick and dirty test for a confident ref/ref site
        # without having access to allele balance info.
        rowSums(as.matrix(tab[,sc.alts]) == 0) > 0 &
        tab$chr == snakemake@wildcards[['chr']],]
write.table(candidate.mosaics[,c(1:2, 4:5)], sep="\t",
    quote=FALSE, row.names=FALSE, file=snakemake@output[[1]])
