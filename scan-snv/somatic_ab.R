args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4)
    stop("usage: balprof.R somatic_positions.txt training.rda fit.rda output.rda\n")

infile <- args[1]
rda <- args[2]
fitfile <- args[3]
outfile <- args[4]

if (file.exists(outfile))
    stop(sprintf("output file %s already exists, please delete it first", outfile))

source("predict.R")

som.tab <- read.table(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE)
str(som.tab)

load(rda)      # loads 'data'
load(fitfile)  # loads 'fits'

somatic.ab.bychrom <- lapply(1:22, function(chrom) {
    hsnps <- data[data$chr==chrom,]
    fit <- as.data.frame(fits[chrom,,drop=FALSE])
    if (ncol(fit)==5)
        colnames(fit) <- c('a', 'b', 'c', 'd', 'logq')

    pos <- som.tab[som.tab$chr == chrom,]$pos
    if (length(pos) == 0) {
        cat(sprintf("no sites to infer on chr%d\n", chrom))
        return(c())
    }

    cat(sprintf("inferring for %d sites on chr%d:%d-%d\n", 
        length(pos), chrom, min(pos), max(pos)))
    system.time(z <- infer.gp(ssnvs=data.frame(pos=pos),
                              fit=fit, hsnps=hsnps, chunk=1, flank=2e5))
    cbind(chr=chrom, pos=pos, z)
})

somatic.ab <- do.call(rbind, somatic.ab.bychrom)

save(somatic.ab, file=outfile)
