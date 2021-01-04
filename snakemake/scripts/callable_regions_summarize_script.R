library(viridis)
sc.sample <- snakemake@wildcards[['sample']]
bulk.sample <- snakemake@config[['bulk_sample']]

# format: 3 columns: chr pos depth (as an integer)
# setting colClass="NULL" prevents reading the column, which saves
# about 2/3rds RAM usage.
sc <- read.table(snakemake@input[['sc']],
    sep='\t', header=TRUE, check.names=FALSE,
    colClasses=c('NULL', 'NULL', 'integer'))
bulk <- read.table(snakemake@input[['bulk']],
    sep='\t', header=TRUE, check.names=FALSE,
    colClasses=c('NULL', 'NULL', 'integer'))

# extra sanity checks: make sure samples names match
if (colnames(sc)[1] != sc.sample)
    stop(sprintf("expected column 3 to be single cell %s, got %s instead",
        sc.sample, colnames(sc)[1]))
if (colnames(bulk)[1] != bulk.sample)
    stop(sprintf("expected column 3 to be bulk sample %s, got %s instead",
        bulk.sample, colnames(bulk)[1]))


# clamp the maximum depth value at 500.
clamp.dp <- 500
x <- data.frame(pmin(clamp.dp, sc[,1]), pmin(clamp.dp, bulk[,1]))
colnames(x) <- c(colnames(sc)[1], colnames(bulk)[1])

# by adding points from (0,0) ... (500,500), the table will be at
# least 500x500 then we can subset to whatever is reasonable.
# subtract the diagnoal after to remove these.
dummy <- data.frame(0:clamp.dp, 0:clamp.dp)
colnames(dummy) <- c(colnames(sc)[1], colnames(bulk)[1])
x <- rbind(x, dummy)

dptab <- table(x)
dptab <- dptab[1:(clamp.dp+1),1:(clamp.dp+1)] - diag(clamp.dp + 1) # diag removes the dummy points

# don't make plots for each chunk; final one is sufficient
#genome.region <- snakemake@config[['gatk_regions']][as.integer(snakemake@wildcards[['gatk_chunk']])]

#pdf(snakemake@output[['pdf']])
#filled.contour(x=0:200, y=0:200, dptab, color.palette=viridis,
    #xlab=sc.sample, ylab=bulk.sample, main=genome.region)
#dev.off()

save(clamp.dp, dptab, file=snakemake@output[[1]])
