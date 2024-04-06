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
            snakemake@input['inttab'],
            snakemake@params['bulk_sample'],
            snakemake@params['sc_sample'],
            snakemake@params['genome'],
            snakemake@params['chrom'],
            snakemake@output['rda'],
            snakemake@params['n_tiles'],
            snakemake@params['n_inits'],
            n_cores=snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 8) {
    cat('if n.cores is unspecified, default is 1\n')
    stop("usage: abmodel_scatter_chrom_tiled_script.R integrated_table.tab.bgz bulk_sample_id sc_sample_id genome chromosome out.rda n.tiles n.inits [n.cores]")
}

int.tab.bgz <- args[1]
bulk.sample <- args[2]
sc.sample <- args[3]
genome <- args[4]
chrom <- args[5]
out.rda <- args[6]
n.tiles <- as.integer(args[7])
n.inits <- as.integer(args[8])

library(future)
if (length(args) == 9) {
    plan(multicore)
    n.cores <- as.integer(args[9])
} else {
    plan(sequential)
    n.cores <- 1
}

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))


suppressMessages(library(scan2))
suppressMessages(library(future))
suppressMessages(library(future.apply))

cat(future::availableCores(), "cores detected by library(future),", n.cores, "requested by user", "\n")
plan(multicore, workers=n.cores)

region <- as(genome.string.to.seqinfo.object(genome)[chrom,], 'GRanges')
hsnps <- read.training.hsnps(path=int.tab.bgz, sample.id=sc.sample, region=region)
# We don't expose hsnp.tilesize to the user.
hsnps <- abmodel.downsample.hsnps(hsnps, hsnp.tilesize=100, n.tiles=n.tiles, verbose=TRUE)

# Context for the C laplace approximation interface
# hsnp.chunksize is the same as hsnp.tilesize (just a poorly chosen name)
ctx <- abmodel.approx.ctx(x=hsnps$pos, y=hsnps$phased.hap1,
                d=hsnps$phased.hap1 + hsnps$phased.hap2,
                hsnp.chunksize=100)

# N.B. we usually show log10(b) and log10(d), but the log likelihood approximator
# expects these in their natural scale.
#
# There are many things that should be done to improve this optimization:
#
#   1. analytical gradient derivation
#
# Algorithm 5.1 in Rasmussen and Williams shows how to compute the partial derivatives
# for an example binary GP. Complexity is O(n^3).
#
#   2. better starting.params
#
# Whether the optimization procedure discovers parameters as good as the grid fit depends
# on the starting parameters. This is, of course, unsettling, but is probably best solved
# by simply running the procedure N times for N sets of starting parameters and choosing
# the best fit.
#
#   3. proper normalization of parameter scale
#
# (3) is really important. R's optim() documentation states that parameters should be
# scaled such that one unit change in a single parameter produces approx. a unit change
# in the log likelihood. This is not even remotely the case for our parameters, especially
# not for the noise parameters (a, c) vs. the distance parameters (b, d). This is currently
# addressed with parameter scales (parscale) and was only roughly worked out by trial and
# error on real data. There are analytical ways to do this better.
#
# Another important point regards the structure of the covariance function: it is a simple
# sum of two identically structured Gaussian kernels, each with 2 parameters. (a,b) are for
# kernel 1 and (c,d) are for kernel 2. Because the two kernels are identical, the parameters
# a and c are not distinguishable (and the same goes for b and d). When grid fitting, we
# artificially choose b <= d by swapping values whenever b>d, but doing this in optim()
# would likely lead to issues with
# taking derivatives. Instead, we supply different parameter scales (parscale) to scale 'd'
# to be ~10-fold larger than b. This is based solely on an a priori expectation that the
# smaller kernel, which usually describes covariance at the scale = read length (~150), is about
# an order of magnitude smaller than the bigger kernel, which usually captures covariance
# at the scale of an amplicon fragment (1-5kb for PTA, longer for MDA).

# these work
#lower.bounds <- c(a=-7,  b=10^2,   c=-7,   d=10^2)
#upper.bounds <- c(a=2,   b=10^5,   c=2,    d=10^5)
#starting.params <- rowMeans(cbind(lower.bounds, upper.bounds))
# experiment
fits <- future.apply::future_lapply(1:n.inits, function(i) {
    lower.bounds <- c(a=-7,  b=10^2,   c=-7,   d=10^2)
    upper.bounds <- c(a=2,   b=10^6,   c=2,    d=10^6)
    #starting.params <- c(0, 150, 0, 1000)
    set.seed(i)
    # choose the random starting params in log-scale for b and d or else
    # the random values will be strongly biased toward large values.
    starting.params <- runif(n=4, min=c(-7,2,-7,2), max=c(2,6,2,6))
    starting.params[2] <- 10^starting.params[2]
    starting.params[4] <- 10^starting.params[4]
    if (starting.params[2] > starting.params[4])
        starting.params[c(2,4)] <- starting.params[c(4,2)]  # b < d

    runtime <- system.time(optimized.params <- optim(par=starting.params,
        fn=function(ps) abmodel.approx.logp(a=ps[1], b=ps[2], c=ps[3], d=ps[4], ctx=ctx, verbose=F),
        method='L-BFGS-B',
        lower=lower.bounds,
        upper=upper.bounds,
        control=list(fnscale=-1, REPORT=1, trace=1, parscale=c(1, 1000, 1, 10000))))
    cat('runtime:\n')
    print(runtime)
    
    pars <- optimized.params$par
    # Enforce b < d
    if (pars[2] > pars[4])
        pars <- pars[c(3,4,1,2)]
    
    # match the format of other param estimator
    fit <- data.frame(iter=i, init.a=starting.params[1], init.b=starting.params[2],
            init.c=starting.params[3], init.d=starting.params[4],
            a=pars[1], b=pars[2], c=pars[3], d=pars[4], logp=optimized.params$value)
    list(fit, optimized.params)
})

save(n.tiles, chrom, fits, file=out.rda, compress=FALSE)

if ('snakemake' %in% ls()) {
    sink()
}
