library(viridis)
sc.sample <- snakemake@wildcards[['sample']]
bulk.sample <- snakemake@config[['bulk_sample']]

# format: 6 columns: chr pos singlecell chr pos bulk
x <- read.table(snakemake@input[[1]], sep='\t', header=TRUE, check.names=FALSE)

# extra sanity checks: make sure columns 3 and 6 are the correct
# single cell and bulk samples
if (colnames(x)[3] != sc.sample)
    stop(sprintf("expected column 3 to be single cell %s, got %s instead",
        sc.sample, colnames(x)[3]))
if (colnames(x)[6] != bulk.sample)
    stop(sprintf("expected column 6 to be bulk sample %s, got %s instead",
        bulk.sample, colnames(x)[6]))


# by adding points from (0,0) ... (200,200), the table will be at
# least 200x200 # then we can subset to whatever is reasonable.
# subtract the diagnoal after to remove these.
dummy <- cbind(0:200,0:200)
colnames(dummy) <- colnames(x)[c(3,6)]
x <- rbind(x[,c(3,6)], dummy)
dptab <- table(x)
dptab <- dptab[1:201,1:201] - diag(201) # diag removes the dummy points

# don't make plots for each chunk; final one is sufficient
#genome.region <- snakemake@config[['gatk_regions']][as.integer(snakemake@wildcards[['gatk_chunk']])]

#pdf(snakemake@output[['pdf']])
#filled.contour(x=0:200, y=0:200, dptab, color.palette=viridis,
    #xlab=sc.sample, ylab=bulk.sample, main=genome.region)
#dev.off()

save(dptab, file=snakemake@output[['rda']])
