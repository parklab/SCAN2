#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) {
        ret <- unlist(c(
            snakemake@input['scan2_object'],
            snakemake@input['spatial_depth'],
            snakemake@input['spatial_abmodel'],
            snakemake@output['rda']
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
    stop("usage: spatial_sens_integrate.R scan2_object.rda spatial_sens_depth.txt spatial_sens_abmodel.txt output.rda")
}

scan2.path <- args[1]
spatial.depth.path <- args[2]
spatial.abmodel.path <- args[3]
out.rda <- args[4]

for (f in out.rda)
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))

suppressMessages(library(scan2))


load(scan2.path) # loads "results"

cat("Reading spatial depth covariates..\n")
depth.covs <- fread(spatial.depth.path)
depth.covs$chr = as.character(depth.covs$chr)
cat("Reading spatial abmodel covariates..\n")
abmodel.covs <- fread(spatial.abmodel.path)
abmodel.covs$chr = as.character(abmodel.covs$chr)

cat("Integrating..\n")
results <- integrate.spatial.sensitivity.covariates(object=results,
    abmodel.covs=abmodel.covs, depth.covs=depth.covs, sens.tilewidth=1e3)

cat("Writing results to", out.rda, "\n")
save(results, file=out.rda, compress=FALSE)

if ('snakemake' %in% ls()) {
    sink()
}
