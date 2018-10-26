args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 1 & length(args) != 2)
    stop("usage: check_grids gridfit_dir [gridn]")

dir=args[1]
gridn=4
if (length(args) == 2)
    gridn = as.integer(args[2])

source("~/balance/gridfit_slurm/bounds.R")
rf <- read.fit(gridn, dir=dir, 20)
