# General routines for calculating enrichment over a set of intervals
# defined by a GenomicRanges object.
# Usually these objects will come from BED files; gene models are the
# primary example of an exception.


use.future <- TRUE

# gmuts is a GRanges-based mutation table
# XXX: It is a ~10-fold speedup to annotate all permuted sites in one
#      call to findOverlaps. There must be some up-front cost to running
#      findOverlaps.
# NOTE NOTE NOTE NOTE ----
#   this requires a disjoint BED so that no permuted mutation can map
#   to multiple features.
# NOTE NOTE NOTE NOTE ----
# --- Update: this is now implemented
map.feature <- function(gmuts, edata, paranoid=FALSE) {
    # 99% of time is spent in findOverlaps
    # for some reason, it is 10x faster to run findOverlaps on a
    # c()ed set of GRanges than it is to call individually on each
    # GRange.
    #ols <- findOverlaps(edata$gbed, gmuts, type='any')
    #mcols(edata$gbed[queryHits(ols)])[[attr(edata$gbed, 'feature.name')]]
    ols <- findOverlaps(gmuts, edata$gbed, type='any')
    # currently paranoid=T is the default, will change after more testing
    if (paranoid) {
        cat('paranoid=TRUE: ensuring all FROM mappings unique:\n')
        print(system.time(anydups <- any(duplicated(from(ols)))))
        if (anydups)
            stop('error: one mutation mapped to multiple BED intervals; is your BED disjoint?')
    }
    # data.table provides some important speed benefits for later grouping
    data.table(from=from(ols), to=features(edata=edata)[to(ols)])
    #mcols(edata$gbed[to(ols)])[[attr(edata$gbed, 'feature.name')]]
}

# Get counts for all possible feature types (including features
# with 0 count).
# x is a vector of the feature being measured (i.e., the return
# value of map.feature()).
count.feature <- function(x, edata) {
    # pre-factor code: O(n^2), can't be used for very large
    # feature sets like genes.
    #sapply(attr(edata$gbed, 'atypes'), function(ty) sum(x == ty))
    ret <- tabulate(x, nbins=nlevels(x))
    names(ret) <- levels(x)
    ret
}

features <- function(en, edata) {
    if (!missing(en))
        edata <- en$edata

    mcols(edata$gbed)[[attr(edata$gbed, 'feature.name')]]
}

# Interface for enrichment analysis. At the minimum, requires
# counting and mapping functions. Additional data (such as a BED
# file of regions) can be specified in '...'.
enrich.data <- function(count.fn=count.feature, map.fn=map.feature, ...) {
    list(count.fn=count.fn, map.fn=map.fn, ...)
}

# Standardize to adding 'chr' tag in sequence names
read.feature.bed <- function(bedfile, feature.name='feature', remove.chr.prefix=FALSE) {
    bed <- read.table(bedfile, stringsAsFactors=F)
    colnames(bed)[1:4] <- c('chr', 'start', 'end', feature.name)
    # Use GRanges for faster intersecting
    # Add 'chr' if missing.
    if (remove.chr.prefix)
        bed$chr <- sub('chr', '', bed$chr)
    ret <- GRanges(seqnames=ifelse(grepl('chr', bed$chr) | remove.chr.prefix,
                bed$chr, paste0('chr', bed$chr)),
        ranges=IRanges(start=bed$start, end=bed$end))
    #if (!isDisjoint(ret))
        #stop('supplied BED must not contain overlapping intervals')

    mcols(ret)[[feature.name]] <- factor(bed[,4])
    # I hate using factors, but this is the only way I know to get
    # tabulate() to always count hits to every possible level and to
    # report it in the same order every time.
    attr(ret, 'feature.name') <- feature.name
    # no longer necessary using factors
    #attr(ret, 'atypes') <-
        #unique(mcols(ret)[[attr(ret, 'feature.name')]])
    ret
}

# For BEDs where the regions indicate an "ON" state: e.g.,
# ChIP histone signals, DNA binding protein peaks, etc.
# In these BEDs, the "OFF" states are the intervals not
# covered by the BED file.
# inout BEDs do not have any additional score metadata.
read.inout.bed <- function(bedfile, feature.name='on', add.chr.prefix=TRUE) {
    cat('WARNING: read.inout.bed assumes UCSC.hg19\n')
    bed <- read.table(bedfile, stringsAsFactors=F)
    if (ncol(bed) == 4)
        colnames(bed)[1:4] <- c('chr', 'start', 'end', feature.name)
    # allow 3 column bed for convenience; assume there is only one state
    if (ncol(bed) == 3)
        colnames(bed)[1:3] <- c('chr', 'start', 'end')
    # Use GRanges for faster intersecting
    # Add 'chr' if missing.
    seqinfo <- seqinfo(getBSgenome(BSgenome.Hsapiens.UCSC.hg19))[paste0('chr',c(1:22,'X','Y'))]
    if (!add.chr.prefix) 
        seqnames(seqinfo) <- sub('chr','',seqnames(seqinfo))

print(seqinfo)

    on <- GRanges(seqnames=ifelse(grepl('chr', bed$chr) | !add.chr.prefix,
                bed$chr, paste0('chr', bed$chr)),
        ranges=IRanges(start=bed$start, end=bed$end),
        seqinfo=seqinfo)

    mcols(on)[[feature.name]] <- factor(x='in', levels=c('in','out'))
    # automatically determining the "off" intervals is what requires
    # the BSgenome's seqinfo.
    # similar to gaps(), but does not do weird things with +/- strand
    off <- setdiff(as(seqinfo(on), 'GRanges'), on)
    mcols(off)[[feature.name]] <- 'out'
    ret <- c(on, off)
    #if (!isDisjoint(ret))
        #stop('supplied BED must not contain overlapping intervals')
    # I hate using factors, but this is the only way I know to get
    # tabulate() to always count hits to every possible level and to
    # report it in the same order every time.
    #attr(ret, 'feature.name') <- feature.name
    attr(ret, 'feature.name') <- feature.name
    #attr(ret, 'atypes') <-
        #unique(mcols(ret)[[attr(ret, 'feature.name')]])
    ret
}

# Just translate a mutation table into a GRanges object
gr <- function(muts, add.chr.prefix=FALSE) {
    g <- GRanges(seqnames=ifelse(grepl('chr', muts$chr) | !add.chr.prefix, 
        muts$chr, paste0('chr', muts$chr)),
        ranges=IRanges(start=muts$pos, width=1))
    g$type.and.ctx <- muts$type.and.ctx
    g
}

map.and.count <- function(combined.set, edata, verbose=TRUE) {
    if (use.future) {
        stop('use.future is no longer supported')
    } else {
        # jointly map all sites to the BED
        st <- system.time(mapping <- edata$map.fn(combined.set, edata))
        # figure out which permutation ID corresponds to each mapping
        # this is the critical requirement for non-disjoint BEDs
        mapping[,perm.id := combined.set$perm.id[mapping$from]]
        if (verbose) { cat('mapping: '); print(st) }
        st <- system.time(
            ret <- mapping[, as.list(edata$count.fn(to)), by=perm.id][,-'perm.id'])
        if (verbose) { cat('tabulating: '); print(st) }
    }

    ret
}

make.perm.matrix <- function(czperms, edata, verbose=TRUE) {
    if (verbose) cat(sprintf("Analyzing %d permutations: [\n", max(czperms$perm.id)))
    ret <- map.and.count(czperms, edata=edata, verbose=verbose)
    if (verbose) cat('] 100%\n')
    ret
}

bootstrap <- function(gmuts, edata, n.boot, verbose=TRUE) {
    if (verbose) cat(sprintf('Bootstrapping %d samples: [\n', n.boot))

    # the combined dataset has to have an ID called 'perm.id' for
    # map.and.count to group the mapping prior to counting.
    cat('sampling', length(gmuts)*n.boot, 'for bootstrapping\n')
    boot.obs <- gmuts[sample(length(gmuts), size=length(gmuts)*n.boot, replace=TRUE),]
    # Rle is a problem for data.table, unfortunately
    #boot.obs$perm.id <- Rle(values=1:n.boot, lengths=length(gmuts))
    boot.obs$perm.id <- rep(1:n.boot, each=length(gmuts))
    ret <- map.and.count(boot.obs, edata=edata, verbose=verbose)
    if (verbose) cat('] 100%\n')
    ret
}

# old way, prior to massive speedup due to combining GRanges
# and feeding into a single findOverlaps
bootstrap.old <- function(gmuts, edata, n.boot, verbose=TRUE) {
    if (verbose) cat(sprintf('Bootstrapping %d samples: [\n', n.boot))

    work <- function(i) {
        if (verbose) { if (i %% 100 == 0) cat(i) else if (i %% 10 == 0) cat('.') }
        make.real.obs(gmuts[sample(length(gmuts), size=length(gmuts), replace=TRUE),],
                      edata)
    }

    if (use.future) {
        stop('use.future is no longer supported')
        boots <- t(future_sapply(1:n.boot, work, future.seed=1234))
    } else {
        boots <- t(sapply(1:n.boot, work))
    }

    if (verbose) cat('] 100%\n')
    list(boots=boots)
}

# OLD STRATEGY: uses a list of GRanges objects. Turns out that
# it is MUCH faster to use a combined object (I call these
# czperms).
# Make a matrix of counts given a list of zipped permutations.
# (i.e., each permutation is the combined set of permutations
# across samples and they have already been converted into
# GRanges.)
make.perm.matrix.old <- function(zperml, edata, verbose=TRUE) {
    if (verbose) cat(sprintf("Analyzing %d permutations: [", length(zperml)))

    work <- function(i) {
        z <- zperml[[i]]
        if (verbose) { if (i %% 100 == 0) cat(i) else if (i %% 10 == 0) cat('.') }
#cat('\nedata\n')
#print(edata$gbed)
#cat('\ngmuts\n')
#print(z)
#exit(1)
        edata$count.fn(edata$map.fn(gmuts=z, edata=edata), edata=edata)
    }

    if (use.future) {
        ret <- future_sapply(1:length(zperml), work, future.seed=1234)
    } else {
        ret <- sapply(1:length(zperml), work)
    }

    if (verbose) cat('] 100%\n')
    t(ret)
}

# Just for convenience
#make.real.obs <- function(gmuts, edata)
    #edata$count.fn(edata$map.fn(gmuts=gmuts, edata=edata), edata=edata)

# Either standardize measurements (as a z-score) or compute the
# enrichment of each observation (including permutations) relative
# to the mean count of permutations.
transform <- function(e, standardize=FALSE) {
    means <- colMeans(e$perm.mat)
    sds <- apply(e$perm.mat, 2, sd)
    if (standardize) {
        bmat <- apply(e$perm.mat, 2, function(col) (col-mean(col)) / sd(col))
        breal <- (e$real.obs - means) / sds
    } else {
        bmat <- apply(e$perm.mat, 2, function(col) col/mean(col))
        breal <- e$real.obs / means
    }
    ret <- list(real=as.vector(breal), pmat=bmat)
    if ('boots' %in% names(e)) {
        if (standardize)
            bs <- sapply(1:length(means), function(i)
                (e$boots[,i] - means[i]) / sds[i])
        else
            bs <- sapply(1:length(means), function(i)
                e$boots[,i] / means[i])
        ret <- c(ret, list(boots=bs))
    }

    ret
}

# Compute p-value of enrichment by determining how often permuted
# counts are further from 1 than the observed count.
pval <- function(e) {
    te <- transform(e)
    lre <- log2(te$real)
    lpm <- log2(te$pmat)
    sapply(1:length(lre), function(i) mean(abs(lpm[,i]) >= abs(lre[i])))
}


# add.outside - given the total number of mutations analyzed,
#     compute the number of mutations outside of any interval.
#     any positive integer will lead to computation.
add.outside <- function(e, n) {
    #cat('WARNING: bootstraps must be recomputed to get proper values, using subtraction method to estimate instead\n')
    #cat('    followup: actually I think subtraction does produce correct estimates.\n')
    e$perm.mat <- cbind(n - rowSums(e$perm.mat), e$perm.mat)
    e$real.obs <- c(n - sum(e$real.obs), e$real.obs)
    if ('boots' %in% names(e)) {
        # XXX: this is actually not correct
        e$boots <- cbind(n - rowSums(e$boots), e$boots)
    }
    e
}


esummary <- function(e, bootstrap.ci=0.95) {
    pv <- pval(e)
    sumdf <- data.frame(
        #state=names(e$real.obs),
        pval=pv,
        padj=p.adjust(pv, method='holm'),
        fdr=round(p.adjust(pv, method='fdr'),4),
        #qval=round(qvalue(pv, lambda = seq(max(0.05,min(pv)), min(0.95,max(pv)), 0.01))$qvalues, 4),
        enr=round(e$real.obs / colMeans(e$perm.mat), 3),
        obs=e$real.obs,
        perm.mean=colMeans(e$perm.mat),
        perm.med=apply(e$perm.mat, 2, median)
    )
    # add bootstrapping stats: number of bootstraps, 95% LB and UB
    if ('boots' %in% names(e)) {
        ci <- bootstrap.ci(e, ci=bootstrap.ci)
        sumdf$n.bootstraps <- nrow(e$boots)
        lb.name <- paste0('boot.', as.character(bootstrap.ci), '.lb')
        ub.name <- paste0('boot.', as.character(bootstrap.ci), '.ub')
        sumdf[[lb.name]] <- ci[1,]
        sumdf[[ub.name]] <- ci[2,]
    }
    if (!is.null(names(e$real.obs)))
        rownames(sumdf) <- names(e$real.obs)
    sumdf
}


bootstrap.ci <- function(boots, ci=0.95) {
    alpha=1-ci
    apply(boots$boots, 2, quantile, probs=c(alpha/2, 1-alpha/2))
}

# Typical enrichment analysis
# By default, the analysis counts the features defined in 'gbed'.
# Other analyses can be performed by supplying different count.fn
# and map.fn functions.
analyze.enrich <- function(muts, zperml, edata, n.boot=1e4) {
    perm.mat <- make.perm.matrix(zperml, edata=edata)
    # permutations are no longer used, so delete them to save memory
    # (they're huge, several GBs).
    cat('removing permutations to conserve RAM\n')
    cat('XXX: fixme: need to remove from caller env, not this one\n')
    rm(zperml, envir=parent.frame(n=1))
    print(gc())
    #real.obs <- make.real.obs(gr(muts, add.chr.prefix=F), edata=edata)
    gmuts <- gr(muts, add.chr.prefix=F)
    gmuts$perm.id <- 1
    real.obs <- unlist(map.and.count(gmuts, edata=edata, verbose=FALSE))
    boots <- bootstrap(gr(muts, add.chr.prefix=F), edata=edata, n.boot=n.boot)
    list(perm.mat=perm.mat, real.obs=real.obs, edata=edata, boots=boots)
}

# hack for plotting things in a different order
# WARNING: reordering the matrices breaks the ordering between
# levels(features) and the matrices and observations.
reorder <- function(en, new.order) {
    list(perm.mat=en$perm.mat[,new.order],
         real.obs=en$real.obs[new.order],
         boots=en$boots[,new.order],
         edata=en$edata)
}

# allow for subsetting an enrichment object to the set of features 'fts'
subset <- function(en, fts) {
    list(
        perm.mat=en$perm.mat[,fts],
        real.obs=en$real.obs[fts],
        boots=en$boots[,fts],
        edata=en$edata
    )
}

# Plotting methods
# 'en' is the result of analyze.enrich()
# type:
#   'b' - plot lines and boxplots
#   'x' - plot boxplots only
#   'l' - plot lines only
# bootstrap=FALSE suppresses bootstrap plotting even if present
plot.enrich <- function(en, bootstrap=0.95, normalize=FALSE, type='b', lcol=1, lwd=2, ltype='b', ...) {
    te <- transform(en, standardize=normalize)
    ylim <- NA

    boxranges <- c(apply(te$pmat, 2, quantile, prob=0.025),
        apply(te$pmat, 2, quantile, prob=0.975))
    if (!('ylim' %in% names(pairlist(...)))) {
        if (type == 'b')
            ylim <- range(c(te$real, boxranges, bootstrap.ci(te, bootstrap)))
        else if (type == 'x')
            ylim <- range(boxranges)
        else if (type == 'l')
            ylim <- range(c(te$real, bootstrap.ci(te, bootstrap)))
    }
    ylab <- ifelse(normalize, 'Mutation enrichment (Z-score)',
        'Observed/Expected')

    bp <- 1:ncol(te$pmat)
    if (type == 'b' | type == 'x') {
        # range reduces whisker length; intent is to remove whiskers
        if ('ylim' %in% names(pairlist(...)))
            boxplot(te$pmat, outline=F, range=0.0000001, border='grey',
                ylab=ylab, ...)
        else 
            boxplot(te$pmat, outline=F, range=0.0000001, border='grey',
                ylim=ylim, ylab=ylab, lty='solid', pars=list(staplewex=0), ...)
        # draw custom whiskers corresponding to 95% CIs on permutations
        arrows(x0=bp, y0=apply(te$pmat, 2, quantile, prob=0.025),
            x1=bp, y1=apply(te$pmat, 2, quantile, prob=0.975),
            angle=90, code=3, length=0.05, col='grey')
    }
    if (type == 'b' | type == 'l') {
        plotf <- if (type == 'b' | c(pairlist(...)[['add']], FALSE)[1] == TRUE) lines else plot
        if ('ylim' %in% names(pairlist(...)))
            plotf(bp, te$real, type=ltype, lwd=lwd, pch=20, col=lcol, ...)
        else
            plotf(bp, te$real, type=ltype, lwd=lwd, pch=20, ylim=ylim, col=lcol, ...)
        if ('boots' %in% names(en) & bootstrap != FALSE) {
            cis <- bootstrap.ci(te, bootstrap)
            arrows(x0=bp, y0=cis[1,], x1=bp, y1=cis[2,],
                angle=90, code=3, length=0.05, col=lcol)
        }
    }
}

# to change the order of bars in the plot, use reorder() on e first
barplot.enrich <- function(e, ...) {
    es <- esummary(e)

    bp <- barplot(100*(es$enr-1),
        names.arg=rownames(es),
        las=3, col='#4c4c4c',
        border=NA, 
        ylab='Enrichment or depletion (%)', ...)
    abline(h=0)

    bs <- bootstrap.ci(transform(e), 0.95)
    arrows(x0=bp, y0=100*(bs[1,]-1), x1=bp, y1=100*(bs[2,]-1),
        angle=90, code=3, length=0.05,lwd=1)
}

##################################################################################
# Methods for making Volcano-style plots out of multiple enrichment analyses
##################################################################################

# load an enrichment object file and do some additional things:
#    1. reorder the features
#    2. if a compute.outside > 0, then add the 'outside' interval
read.e <- function(file, compute.outside=-1, check.features=NULL) {
    print(file)
    load(file)
    if (!is.null(check.features)) {
        if (!all(unique(features(e)) %in% check.features))
            stop(paste('features in file', file, 'were not in the specified feature set'))
    }
    e <- reorder(e, order(unique(features(e))))
    if (compute.outside > 0) {
        e <- add.outside(e, compute.outside)
    }
    e
}


# just reads and reorders the enrichment objects and ensures
# that all enrichment objects are compatible (e.g., were run
# on the same feature sets).
read.es <- function(files, compute.outside=-1) {
    require(GenomicRanges) # the features() function accesses the GRanges object

    # load just the first file to get the feature set.
    load(files[1])
    feature.set <- unique(features(e))
    es <- lapply(files, function(fn) {
        read.e(fn, compute.outside=compute.outside, check.features=feature.set)
    })
}


# map each enrichment object to a single (enrichment, p-value)
# pair. this requires choosing a specific feature on which to
# compute enrichment.
# produces a matrix suitable for evolcano()
# min.pval - should be 1 / (# of permutations used). Prevents
#    Inf p-values when none of the permutations exceeded the
#    observed values.
volcanoize <- function(es, feature='in peak', min.pval=1e-4) {
    # doing feature by name
    if (is.character(feature)) {
        if (!(feature %in% unique(features(es[[1]]))))
            stop(paste0("feature '", feature, "' not in first enrichment object"))
    } else if (is.integer(feature)) {
        if (feature > length(unique(features(es[[1]]))))
            stop(paste0("feature index ", feature, " exceeds the length of the first enrichment object"))
    }

    ret <- sapply(es, function(e) {
        ret <- unlist(esummary(e)[feature, c('enr', 'pval')])
        ret[2] <- max(min.pval, ret[2])
        ret
    })
    ret[2,] <- -log10(ret[2,] + min.pval)
    colnames(ret) <- names(es)
    rownames(ret) <- c('enr','sig')
    t(ret)
}


evolcano <- function(emat, labels=rep('', nrow(emat)), ...) {
    require(basicPlotteR)

    z <- emat
    # points will be red if they're significant OR labeled
    plot(z[,'enr'], z[,'sig'], pch=16,
        col=ifelse(z[,'sig'] <= 1 & labels == '', 'grey', 1 + (labels != '')),
            ylim=c(0,4.5), xlab='Enrichment (or depletion) ratio (obs/exp)',
            ylab='Significance: -log10(p-value)', cex=2, ...)

    if (any(labels != '')) {
        z <- z[labels != '',]
        addTextLabels(xCoords=z[,'enr'], yCoords=z[,'sig'],
                    labels=labels[labels != ''], col.line=2,
                    col.label=2, cex.label=0.9)
    }
    abline(h=1, lty='dotted')
    abline(v=1, lty='dotted')
}


evolcano.old.encode <- function(df, meta, labels=c('brain', 'all', 'none'), ...) {
    require(basicPlotteR)
    if (ncol(df) != nrow(meta))
        stop('columns of df must correspond exactly to rows of meta')

    z <- data.frame(enr=df[1,], pval=pmin(-log10(min.pval),df[2,]),
        tissue=meta$ANATOMY, specific=meta$STD_NAME,
        stringsAsFactors=FALSE)

    plot(z$enr, z$pval, pch=16,
        col=ifelse(z$pval <= 1 & z$tissue != 'BRAIN', 'grey', 1 + (z$tissue=='BRAIN')),
            ylim=c(0,4.5), xlab='Enrichment (or depletion) ratio (obs/exp)',
            ylab='Significance: -log10(p-value)', ...)

    # keep the parts of the dataframe that are either more significant than 0.1
    # or come from tissues annotated as BRAIN type.
    if (labels != 'none') {
        if (labels == 'all') {
            z <- z[!is.na(z$pval) & (z$pval > 1 | z$tissue=='BRAIN'),]
        } else if (labels == 'brain') {
            z <- z[!is.na(z$pval) & z$tissue=='BRAIN',]
        }
        addTextLabels(xCoords=z$enr, yCoords=z$pval,
                    labels=z$specific, col.line='black',
                    col.label=1 + (z$tissue=='BRAIN'), cex.label=0.9)
    }
    abline(h=1, lty='dotted')
    abline(v=1, lty='dotted')
}



##################################################################################
# For making boxplots using bootstraps
##################################################################################
e.to.bootl <- function(es, feature.idx=2, labels=names(es)) {
    cat("WARNING! WARNING! must modify esummary to not reorder data\n")
    ret <- list(
        enrs=lapply(es, function(e) e$boots[,feature.idx]/mean(e$perm[,feature.idx])),
        pvals=pmax(1e-4,sapply(es, function(e) esummary(e)[feature.idx,]$pval)),
        fdrs=pmax(1e-4, sapply(es, function(e) esummary(e)[feature.idx,]$fdr)),
        labels=labels
        #labels=sub('_promoters','',sub('_enhancers','',labels))
    )
    ret
}

# merge bootls by:
#   1. using median enrichment values as the points for the boxplot and
#   2. using harmonic mean p-values to combine the p-values (for asterisks)
# CAUTION! CAUTION! p-values and FDRs are no longer different after applying
# this function. DO NOT interpret the 'fdrs' output of this function as FDRs.
merge.bootl <- function(bootl, group.labels) {
    harmonicmeanpvs <- sapply(split(bootl$pvals, group.labels),
        function(ps) p.hmp(ps, L=length(ps)))
    list(
        enrs=lapply(split(bootl$enrs, group.labels), function(subenrs) sapply(subenrs, median)),
        pvals=harmonicmeanpvs,
        fdrs=harmonicmeanpvs,
        labels=group.labels[!duplicated(group.labels)])
        
}

boxp <- function(l, add.points=TRUE, asts.at=NA, ...) {
    bp <- boxplot(l$enrs, names=l$labels,
        outline=F, lty='solid', pars=list(staplewex=0),
        las=3, ...)
    abline(h=1, lty='dotted')
    if (!is.na(asts.at)) {
        asts <- sapply(floor(-log10(l$fdrs)),
            function(i) sprintf("%s", paste0(rep('*', i), collapse='')))
        print(asts)
        #text(x=1:length(l$enrs), y=asts.at, labels=asts)
        text(x=1:length(l$enrs), y=bp$stats[5,], labels=asts, adj=c(0.5,0))
    }
    if (add.points) {
        stripchart(l$enrs, vertical=T, pch=20, add=TRUE, method='jitter')
    }
}



# A standard command line analysis:
#   Given a .rda file containing mutations (arg 2) and a corresponding
#   set of random permutations (arg 3), determine the enrichment of
#   mutations over some set of genomic regions defined by the object
#   returned by 'init.enrich'. Results are written to the output .rda
#   file in (arg 4).
#
#   Significance is also assessed by bootstrapping with `n.boot`
#   resamplings.
#
#   (arg 1) specifies the number of cores to use for permutation and
#   bootstrapping analysis.
#
# 'init.enrich' is a function that takes an arbitrary number of
# command line arguments (starting with arg 5) and returns an
# enrichment data object.
# args - allow the caller to override args. this allows the caller to
#        implement different command line options and then reshape the
#        args in the way command.line.analysis expects.
command.line.analysis <- function(init.enrich, args=commandArgs(trailingOnly=TRUE)) {
    #args <- commandArgs(trailingOnly=TRUE)
    if (length(args) < 5) {
        cat("Got command args:\n")
        print(args)
        stop("usage: [enrichment_analysis.r] n.bootstraps n.cores varname1:muts.rda varname2:perms.rda[:use_N] outfile.rda [ optional: input1 [ input2 ... ] ]")
    }

    n.boot <- as.integer(args[1])
    n.cores <- as.integer(args[2])
    mutarg <- args[3]
    permarg <- args[4]
    outfile <- args[5]
    inputdata <- args[6:length(args)]

    mut.varname <- strsplit(mutarg, ":")[[1]][1]
    mutfile <- strsplit(mutarg, ":")[[1]][2]
    perm.params <- strsplit(permarg, ":")[[1]]
    perm.varname <- perm.params[1]
    permfile <- perm.params[2]
    if (length(perm.params) == 3) {
        cat('WARNING: useN: this can fail if a large number of permutations is requested\n')
        perms.useN <- as.integer(perm.params[3])
    } else
        perms.useN <- Inf

    if (file.exists(outfile))
        stop(sprintf("output file %s already exists, please delete it first",
            outfile))

    # do basic arg checking first - loading these libraries takes a while
    suppressMessages(library(GenomicRanges))
    suppressMessages(library(parallel))
    suppressMessages(library(future.apply))
    suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    suppressMessages(library(data.table))

    eobject <- init.enrich(inputdata)

    cat(sprintf("Using %d cores.\n", n.cores))
    # allow objects up to ~10GB; 10,000 permutation passA zperml with mut ctx is 7.70 GB
    # XXX: sometimes future significantly overestimates the size of objects.
    options(future.globals.maxSize=10e9)  
    if (n.cores > 1) {
        cat('    -> using futures for parallelization\n')
        plan(multicore, workers=n.cores)
    } else {
        plan(sequential)
        use.future <<- FALSE
    }

    load(mutfile, verb=T)
    muts <- get(mut.varname)
    load(permfile, verb=T)
    zperml <- get(perm.varname)
    
    
    actual.useN <- max(min(length(zperml), perms.useN), 1)
    cat(sprintf("Using the first %d permutations out of %d total (useN=%g)\n", 
        actual.useN, length(zperml), perms.useN))
    size.of.one.perm <- sum(zperml$perm.id==1)  # all permutations must be the same size
    if (is.finite(perms.useN))
        zperml <- zperml[1:(size.of.one.perm*actual.useN)]
    # else don't do anything, because creating a large index vector can
    # cause a memory error. 

    # 'chr' prefixes are becoming such an issue that we need to check
    # before waiting an hour for this to run.
    cat("chr prefix in mutations?", muts$chr[1], '\n')
    cat("chr prefix in permutations?\n")
#str(seqnames(zperml[[1]]))
str(seqnames(zperml))
    cat("chr prefix in bed file?\n")
str(seqnames(eobject$gbed))

    # hack to remove the other zipped list. these objects are
    # large (~1GB per list) and have to be copied to all fork()ed
    # workers.
    # XXX: not sure why, but this doesn't actually reduce mem usage
    #to.rm <- ifelse(perm.varname == 'nzperml', 'gzperml', 'nzperml')
    #cat(sprintf('HACK: deleting %s to save memory\n', to.rm))
    #rm(to.rm)
    cat('Memory profile before running analysis:\n')
    print(gc())
    e <- analyze.enrich(muts, zperml, eobject, n.boot=n.boot)
    cat('Memory profile after running analysis:\n')
    print(gc())

    save(e, file=outfile)
}

