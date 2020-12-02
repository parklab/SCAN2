data <- read.table(snakemake@input[[1]],
    header=TRUE, stringsAsFactors=FALSE,
    colClasses=c(chr='character'))
save(data, file=snakemake@output[[1]])
