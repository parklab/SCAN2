# read a grid output by the two component model (4 parameters)
# in this case, the parameter values are generated randomly by the
# C program and the output binary is structured differently.
read.grid4 <- function(fname) {
    f <- file(fname, 'rb')
    gn <- readBin(f, integer(), 1)
    m <- t(matrix(readBin(f, numeric(), 5*gn), nrow=5))
    colnames(m) <- c('a', 'b', 'c', 'd', 'logq')
    close(f)
    m
}

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


args <- commandArgs(TRUE)
if (length(args) != 3 & length(args) != 4) {
    stop("Rscript combine.R grid.dir gridn ngrids [output.rda]")
}

grid.dir <- args[1]
gridn <- as.integer(args[2])
ngrids <- as.integer(args[3])

write.to <- NA
if (length(args) == 4)
    write.to <- args[4]

x <- read.fit.one(gridn=gridn, dir=grid.dir, n=ngrids)

if (is.na(write.to)) {
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
} else {
    fits <- best.fit(list(x))
    save(fits, file=write.to)
}
