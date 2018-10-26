# R functions for reading and writing binary file formats
# created by the C mkl-gridfit programs.
read.grid <- function(fname, two.component=FALSE) {
    f <- file(fname, 'rb')
    gn <- readBin(f, integer(), 1)
    a.grid <- readBin(f, numeric(), gn)
    b.grid <- readBin(f, numeric(), gn)
    sn <- readBin(f, integer(), 1)
    a.subset <- readBin(f, integer(), sn)+1  # 0-based indexing
    b.subset <- readBin(f, integer(), sn)+1
    v <- readBin(f, numeric(), sn^2)
    m <- t(matrix(v, ncol=sn))
    close(f)
    list(a.grid=a.grid, b.grid=b.grid, a.subset=a.subset, b.subset=b.subset, Lgrid=m)
}

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

# g is the object returned by read.bin above, specifically in
# two.component=TRUE mode.
melt.grid <- function(g) {
    stop("no longer used")
    data.frame(a=rep(g$a.grid, each=length(g$b.grid)*length(g$c.grid)),
        b=rep(rep(g$b.grid, each=length(g$c.grid)), times=length(g$a.grid)),
        c=rep(g$c.grid, times=length(g$a.grid)*length(g$b.grid)),
        d=g$d,
        L=as.vector(g$Lgrid))
}

# write out the position, hap1, dp info for consumption by the C
# mkl-gridfit program.  format:
# N <integer> number of rows
# pos N*<integer> vector of positions
# hap1 N*<integer> vector of hap1s
# dp N*<integer> vector of DPs
write.bin <- function(df, fname) {
    f <- file(fname, 'wb')
    writeBin(nrow(df), f)
    writeBin(df$pos, f)
    writeBin(df$hap1, f)
    writeBin(df$hap1 + df$hap2, f)
    close(f)
}

read.bin <- function(fname) {
    f <- file(fname, 'rb')
    n <- readBin(f, integer(), 1)
    pos <- readBin(f, integer(), n)
    hap1 <- readBin(f, integer(), n)
    dp <- readBin(f, integer(), n)
    close(f)
    data.frame(pos=pos, hap1=hap1, hap2=dp-hap1)
}
