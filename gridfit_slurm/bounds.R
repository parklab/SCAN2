source("~/balance/utils/binio.R")


# the covariance functions i chose were usually non identifiable,
# to make them identifiable we force the parameters b < d.
read.fit <- function(gridn, dir=".", n=25, chrs=22) {
    lapply(1:chrs, function(chr) read.fit.one(gridn=gridn, dir=sprintf("%s/chr%d", dir, chr), n=n))
}

read.fit.one <- function(gridn, dir, n) {
        x <- do.call(rbind, lapply(1:n, function(g)
            tryCatch({
                fname <- sprintf("%s/g%d.%d.bin", dir, gridn, g)
                #cat(sprintf('reading %s\n', fname))
                read.grid4(fname)
            },
            error = function(e) {
                cat(sprintf("binary file %s could not be read\n", fname))
        })))

        # swapping columns (1,2) and (3,4) to force b < d.
        # ifelse returns a value the same shape as the first argument
        logi.mat <- matrix(rep(x[,2] < x[,4], times=5), ncol=5)
        x <- ifelse(logi.mat, x, x[,c(3,4,1,2,5)])
        x[order(x[,5], decreasing=TRUE),]
}

best.fit <- function(rf) {
    t(sapply(rf, function(x) x[1,]))
}


build.bounds <- function(rf, top.n=50, chrs=1:22) {
    b <- cbind(chrs,
        t(sapply(rf, function(y) apply(y[1:min(nrow(y), top.n),-5], 2, range)))
    )
    b[,4:5] <- log10(b[,4:5])
    b[,8:9] <- log10(b[,8:9])
    b
}

# input is the return value of read.fit()
write.bounds <- function(rf, file, top.n=50, chrs=1:22) {
    b <- build.bounds(rf, top.n=top.n, chrs=chrs)
    write.table(b, file=file, quote=FALSE, row.names=FALSE, col.names=FALSE)
}
