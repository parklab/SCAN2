# updated from reading output of R pipeline to reading output
# of C pipeline:
#   - input files are now bin files, read by binio.R
#   - now uses the 4 parameter model with gaussian kernels
#   - since in 4-d space, contour plots aren't that useful anymore
source("~/balance/utils/binio.R")          # need read.grid4
source("~/balance/gridfit_slurm/bounds.R")
args <- commandArgs(TRUE)
if (length(args) != 3) {
    stop("Rscript combine.R grid.dir gridn ngrids")
}

grid.dir <- args[1]
gridn <- as.integer(args[2])
ngrids <- as.integer(args[3])


x <- read.fit.one(gridn=gridn, dir=grid.dir, n=ngrids)

# the chromosome is ignored here anyway
b <- build.bounds(list(x), chrs=0)

# output is: one new boundary per line
# a pair of lines defines the upper and lower bounds for a parameter
cat(sprintf("%0.10f\n", b[2]))
cat(sprintf("%0.10f\n", b[3]))
cat(sprintf("%0.10f\n", b[4]))
cat(sprintf("%0.10f\n", b[5]))
cat(sprintf("%0.10f\n", b[6]))
cat(sprintf("%0.10f\n", b[7]))
cat(sprintf("%0.10f\n", b[8]))
cat(sprintf("%0.10f\n", b[9]))
