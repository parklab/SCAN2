# Allocate all necessary vectors for the C fitting code.
# Trying to avoid memory copying since these matrices can be large.
abmodel.approx.ctx <- function(x, y, d, hsnp.chunksize=100) {
    n <- hsnp.chunksize

    # numeric() vectors are initialized to 0 by default.
    list(hsnp.chunksize=as.integer(hsnp.chunksize),
         x=x,
         y=y,
         d=d,
         U=numeric(n),
         V=numeric(n),
         B=numeric(n),
         sqrtW=numeric(n),
         K=numeric(n*n),
         A=numeric(n*n))
}

# ctx is the set of working memory space allocated above
abmodel.approx.logp <- function(a, b, c, d, ctx,
    max.it=as.integer(50), verbose=FALSE) {
    result <- .Call("laplace_approx_chunk_cpu",
        ctx$hsnp.chunksize, c(a, b, c, d),
        ctx$x, ctx$y, ctx$d, as.integer(length(ctx$x)),
        ctx$U, ctx$V, ctx$B, ctx$sqrtW, ctx$K, ctx$A,
        max.it, verbose,
        PACKAGE="scansnv")
    return(result)
}

abmodel.sample <- function(n=1000, alim=c(-7,2), blim=c(2,4), clim=c(-7,2), dlim=c(2,6),
    ctx, seed=0, max.it=50, verbose=FALSE) {

    set.seed(seed)
    max.it <- as.integer(max.it)

    params <- data.frame(a=runif(n, min=alim[1], max=alim[2]),
        b=10^runif(n, min=blim[1], max=blim[2]),
        c=runif(n, min=clim[1], max=clim[2]),
        d=10^runif(n, min=dlim[1], max=dlim[2]))

    logps <- mapply(abmodel.approx.logp, a=params$a, b=params$b, c=params$c, d=params$d,
        MoreArgs=list(ctx=ctx, max.it=max.it, verbose=verbose))

    return(cbind(params, logp=logps))
}


K.func <- function(x, y, a, b, c, d) exp(a - (x - y)^2 / b^2) + exp(c - (x-y)^2 / d^2)

# From Rasmussen & Williams 2006.  Calculates the conditional distn of
# the latent GP given the observed points without inverting K.  Since
# the latent GP is MVN, "computing the distribution" only requires
# solving for the mean and covariance.
alg3.2.2 <- function(a, b, c, d, ctx, Xnew) {
    # (Using my notation): this conditional distribution is B|Y,X,D.
    # Approximated by the Laplace method, B|Y,X,D ~ MVN(mode, covmat).
    # mode is the maximizer of the nonapproximate distn and
    # covmat=(K + W^-1)^-1, where W is the Hessian of log p(Y|B).
    mode <- ctx$B[1:length(ctx$x)]
    # this is stored in ctx, but I'm not sure that the C code creates
    # a matrix in the format R expects.
    K <- outer(ctx$x, ctx$x, K.func, a=a, b=b, c=c, d=d)
    covK <- outer(ctx$x, Xnew, K.func, a=a, b=b, c=c, d=d)

    # Infer the GP at some new positions Xnew
    # ( Y - d...) is del log p(Y|B=mode)
    mean.new <- t(covK) %*% (ctx$y - ctx$d * exp(mode) / (1 + exp(mode)))

    # v satisfies: v^T v = k(X_sSNV, X)^T (K + W^-1)^-1 k(X_sSNV, X)
    W <- ctx$d * exp(mode) / (1 + exp(mode))^2
    sqrtW <- sqrt(W)
    L <- t(chol(diag(length(ctx$x)) + outer(sqrtW, sqrtW) * K))
    v <- forwardsolve(L, outer(sqrtW, rep(1, ncol(covK))) * covK)

    cov.new <- outer(Xnew, Xnew, K.func, a=a, b=b, c=c, d=d) - t(v) %*% v

    list(mean=mean.new, cov=cov.new)
}

# form a large enough block around the set of variants "vars",
# infer the Laplace-approximate distribution of B|Y at the training
# sites within the block, then infer the same approximate distribution
# on B*|Y*, the balances at the candidate variant sites.
# returns the mean and variance of the GP at the candidate sites.
infer.gp.block <- function(ssnvs, fit, hsnps, ctx, flank=1e5, max.hsnps=150, verbose=FALSE) {
    a <- fit$a
    b <- fit$b
    c <- fit$c
    dparam <- fit$d

    # row index of the (lower, upper) bounds in hsnps, containing
    # ssnvs$pos +/- flank
    left <- findInterval(range(ssnvs$pos)[1], hsnps$pos)
    up <- findInterval(range(ssnvs$pos)[1] - flank, hsnps$pos)
    up <- max(up, left - max.hsnps)
    right <- findInterval(range(ssnvs$pos)[2], hsnps$pos)
    down <- findInterval(range(ssnvs$pos)[2] + flank, hsnps$pos)
    down <- min(down, right + max.hsnps)
    window <- c(up, down)

    d <- hsnps[max(window[1], 1):min(window[2], nrow(hsnps)),]
    if (verbose) {
        print(window)
        cat(sprintf("infer.gp.block: %d nearby hets\n", nrow(d)))
    }

    # approx. distn of B|Y at the training sites
    #z <- alg3.1(a=a, b=b, c=c, d=dparam, X=d$pos, Y=d$hap1, D=d$hap1 + d$hap2)
    ctx$x <- d$pos
    ctx$y <- d$hap1
    ctx$d <- d$hap1 + d$hap2
    abmodel.approx.logp(a=a, b=b, c=c, d=dparam, ctx=ctx)

    # insert the position of the variant to be tested
    z2 <- alg3.2.2(a=a, b=b, c=c, d=dparam, ctx=ctx, Xnew=ssnvs$pos)
    data.frame(gp.mu=z2$mean, gp.sd=sqrt(diag(z2$cov)))
}


# "chunks" here are NOT the 250 hSNP blocks used in parameter fitting.
# "ssnvs" are the candidate sSNVs. the data frame need only have a 'pos'
#        column, but should only contain candidates from one chromosome
# "hsnps" should be the phased hSNPs used for fitting, but again only
#        from one chromosome corresponding to ssnvs.
infer.gp <- function(ssnvs, fit, hsnps, chunk=2500, flank=1e5, max.hsnps=150) {
    nchunks <- ceiling(nrow(ssnvs)/chunk)
    ctx <- abmodel.approx.ctx(c(), c(), c(), hsnp.chunksize=2*max.hsnps + 10)
    do.call(rbind, lapply(1:nchunks, function(i) {
            cat(sprintf("block %d\n", i))
            start <- 1 + (i-1)*chunk
            stop <- min(i*chunk, nrow(ssnvs))
            infer.gp.block(ssnvs[start:stop,,drop=FALSE],
                fit, hsnps, ctx=ctx, flank=flank, max.hsnps=max.hsnps)
        })
    )
}

# predicting each site independently
infer.gp1 <- function(ssnvs, fit, hsnps, flank=1e5, max.hsnps=150) {
    ctx <- abmodel.approx.ctx(c(), c(), c(), hsnp.chunksize=2*max.hsnps + 10)
    do.call(rbind, lapply(1:nrow(ssnvs), function(i) {
            infer.gp.block2(ssnvs[i,,drop=FALSE],
                fit, hsnps, ctx=ctx, flank=flank, max.hsnps=max.hsnps)
    }))
}
