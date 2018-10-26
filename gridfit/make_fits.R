args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 2)
    stop("usage: Rscript make_fits.R dir output.rda")

source("~/balance/gridfit_slurm/bounds.R")

dir=args[1]
outfile=args[2]

if (file.exists(outfile))
    stop(sprintf("output file %s already exists", outfile))

rf <- read.fit(4, dir=dir, n=20)
fits <- best.fit(rf)
save(fits, file=outfile)
