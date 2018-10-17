# These routines implement the same algorithm 3.1 implemented
# by the C mkl-gridfit-gauss program and algorithm 3.2.
# The C program is necessary to perform the grid search for MLEs
# of the covariance function parameters a, b, c and d because so
# many alg3.1 approximations must be run.
# These R functions are only used to predict AB at sSNV candidate
# loci, meaning they do not need to be as optimized.

zzz <- NULL  # for gradient maxing; wish there was a better way
zzzchunk <- list()

K.func <- function(x, y, a, b, c, d) exp(a - (x - y)^2 / b^2) + exp(c - (x-y)^2 / d^2)


# 2017-2-3: also moved from absolute convergence criteria to relative
# 2017-2-3: optimized a few bits of code
#    1. for a diagonal matrix D, DMD = DD^T * M, where * represents
#       elementwise mutliplication.  This is O(n^2), compared to the
#       O(n^3) of matrix multiplication.
#    2. same: MD for matrix M, diagonal matrix D can avoid O(n^3)
#       matrix multiplication
alg3.1 <- function(a, b, c, d, X, Y, D, max.it=50, verbose=TRUE) {
    if (verbose)
        cat(sprintf("alg3.1: (a=%0.5f, b=%0.5f, c=%0.5f, d=%0.5f)\t",
            a, b, c, d))

    n <- length(X)   # will be convenient
    # K.func taken from global context
    K <- outer(X, X, K.func, a=a, b=b, c=c, d=d)
    B <- rep(0, n)
    iter <- 0
    logqs <- -Inf
    repeat {
        iter <- iter + 1
        # The component of the Hessian from log p(Y|B)
        # notice: W is a diagonal matrix
        W <- D * exp(B) / (1 + exp(B))^2
        sqrtW <- sqrt(W)
        L <- t(chol(diag(n) + outer(sqrtW, sqrtW) * K))
        # The term in () is the gradient of log p(Y|B).  In my derivation,
        # grad u = grad{log p(Y|B)} - K^-1 B.  However, since
        # hess u = -W - K^-1, we can see that K^-1 B = hess(u) B + WB, which
        # is where this term arises.
        var1 <- W*B + (Y - D * exp(B) / (1 + exp(B)))
        var2 <- var1 -
          sqrtW * backsolve(t(L), forwardsolve(L, (outer(sqrtW, rep(1, n)) * K) %*% var1))
        var2 <- drop(var2)
        B <- drop(K %*% var2)

        logqs <- c(logqs, drop(-var2 %*% B/2 + (B %*% Y - sum(D * log(1 + exp(B)))) - sum(log(diag(L)))))
        if (abs((logqs[iter+1] - logqs[iter])/logqs[iter+1]) < sqrt(.Machine$double.eps))
            break
        if (iter > max.it)
            stop(sprintf("reached max.it=%d iterations without converging", max.it))
    }
    if (verbose)
        cat(sprintf("\t%d iterations\n", iter))

    # K, K.inv.f (known as a in R&W '06) needed by algorithm 5.1
    retval <- list(iterations=iter, max.it=max.it, B.=B,
                    logq=logqs[iter+1], logqs=logqs, K=K, K.inv.f=var2,
                    a=a, b=b, c=c, d=d)
    zzz <<- retval   # wish there were a better way
    retval
}

# From Rasmussen & Williams 2006.  Calculates the conditional distn of
# the latent GP given the observed points without inverting K.  Since
# the latent GP is MVN, "computing the distribution" only requires
# solving for the mean and covariance.
alg3.2 <- function(a, b, c, d, Xfit, Xnew, Yfit, Dfit, alg3.1ret, Dnew) {
    # (Using my notation): this conditional distribution is B|Y,X,D.
    # Approximated by the Laplace method, B|Y,X,D ~ MVN(mode, covmat).
    # mode is the maximizer of the nonapproximate distn and
    # covmat=(K + W^-1)^-1, where W is the Hessian of log p(Y|B).
    mode <- alg3.1ret$B.
    K <- outer(Xfit, Xfit, K.func, a=a, b=b, c=c, d=d)
    covK <- outer(Xfit, Xnew, K.func, a=a, b=b, c=c, d=d)

    # Infer the GP at some new positions Xnew
    # ( Y - d...) is del log p(Y|B=mode)
    mean.new <- t(covK) %*% (Yfit - Dfit * exp(mode) / (1 + exp(mode)))

    # v satisfies: v^T v = k(X_sSNV, X)^T (K + W^-1)^-1 k(X_sSNV, X)
    W <- Dfit * exp(mode) / (1 + exp(mode))^2
    sqrtW <- sqrt(W)
    L <- t(chol(diag(length(Xfit)) + outer(sqrtW, sqrtW) * K))
    v <- forwardsolve(L, outer(sqrtW, rep(1, ncol(covK))) * covK)

    cov.new <- outer(Xnew, Xnew, K.func, a=a, b=b, c=c, d=d) - t(v) %*% v

    retval <- list(mean=mean.new, cov=cov.new)
    retval
}

# form a large enough block around the set of variants "vars",
# infer the Laplace-approximate distribution of B|Y at the training
# sites within the block, then infer the same approximate distribution
# on B*|Y*, the balances at the candidate variant sites.
# returns the mean and variance of the GP at the candidate sites.
infer.gp.block <- function(ssnvs, fit, hsnps, flank=1e5) {
    a <- fit$a
    b <- fit$b
    c <- fit$c
    dparam <- fit$d

    # row index of the (lower, upper) bounds in hsnps, containing
    # ssnvs$pos +/- flank
    print(range(ssnvs$pos))
    window <- findInterval(range(ssnvs$pos) + c(-flank, flank), hsnps$pos)
    print(window)
    d <- hsnps[max(window[1], 1):min(window[2], nrow(hsnps)),]
    cat(sprintf("infer.gp.block: %d nearby hets\n", nrow(d)))

    # approx. distn of B|Y at the training sites
    z <- alg3.1(a=a, b=b, c=c, d=dparam, X=d$pos, Y=d$hap1, D=d$hap1 + d$hap2)

    # insert the position of the variant to be tested
    z2 <- alg3.2(a=a, b=b, c=c, d=dparam,
        Xfit=d$pos, Yfit=d$hap1, Dfit=d$hap1 + d$hap2,
        Xnew=ssnvs$pos, alg3.1ret=z)
    return(data.frame(gp.mu=z2$mean, gp.sd=sqrt(diag(z2$cov))))
}

# "chunks" here are NOT the 250 hSNP blocks used in parameter fitting.
# "ssnvs" are the candidate sSNVs. the data frame need only have a 'pos'
#        column, but should only contain candidates from one chromosome
# "hsnps" should be the phased hSNPs used for fitting, but again only
#        from one chromosome corresponding to ssnvs.
infer.gp <- function(ssnvs, fit, hsnps, chunk=2500, flank=1e5) {
    nchunks <- ceiling(nrow(ssnvs)/chunk)
    do.call(rbind, lapply(1:nchunks, function(i) {
            cat(sprintf("block %d\n", i))
            start <- 1 + (i-1)*chunk
            stop <- min(i*chunk, nrow(ssnvs))
            infer.gp.block(ssnvs[start:stop,,drop=FALSE],
                fit, hsnps, flank=flank)
        })
    )
}
