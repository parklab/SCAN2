# return strings appropriate for analyzing multiset overlaps with upset
# gatk is a gatk-like matrix with only (genotype, ref, alt) triplets
cell.inclusion <- function(gatk, min.alt=1) {
    cn <- colnames(gatk[,seq(1, ncol(gatk), 3)])
    apply(gatk[,seq(3, ncol(gatk), 3)], 1,
        function(row) paste(cn[!is.na(row) & row >= min.alt], collapse='&'))
}

# given two data frames of somatic and germline locations, annotate
# the somatic data frame with the position of the nearest germline entry.
find.nearest.germline <- function(som, germ) {
    som$nearest.het <- NA
    for (chr in 1:22) {
        gpos <- germ$pos[germ$chr == chr]
        spos <- som$pos[som$chr == chr]
        gidx <- findInterval(spos, gpos)
        nearest.idx <- ifelse(abs(gpos[gidx] - spos) <= abs(gpos[gidx+1] - spos), gidx, gidx+1)
        som$nearest.het[som$chr == chr] <- gpos[nearest.idx]
    }
    som
}

# Reject somatic sites with data points lying outside of the specified
# quantiles (min.q, max.q) at germline sites.
# These types of tests are controlled by SENSITIVITY to germline sites
# (which hopefully applies to somatic sites) rather than the AB tests
# that aim to control the FDR.
test.against.germline <- function(germline, somatic, min.q=0, max.q=1) {
    x <- quantile(germline[!is.na(germline)], probs=c(min.q, max.q))
    min.cutoff <- x[1]
    max.cutoff <- x[2]
    somatic >= min.cutoff & somatic <= max.cutoff
}

read.gatk.tab <- function(filename) {
    read.table(filename, header=T, stringsAsFactors=F)
}

# if n is specified, a random subsample of size n will be written instead
# used for estimating distributions for germline sites
write.gatk.positions <- function(filename, gatk, n=nrow(gatk)) {
    gatk <- gatk[sort(sample(nrow(gatk), size=n)),]
    write.table(cbind(sprintf("%d:%d-%d", gatk$chr, gatk$pos, gatk$pos), gatk$altnt),
        file=filename, quote=F, row.names=F, col.names=F, sep='\t')
}

muttype.map <- c(
    'A>C'='T>G',
    'A>G'='T>C',
    'A>T'='T>A',
    'C>A'='C>A',
    'C>G'='C>G',
    'C>T'='C>T',
    'G>A'='C>T',
    'G>C'='C>G',
    'G>T'='C>A',
    'T>A'='T>A',
    'T>C'='T>C',
    'T>G'='T>G'
)

# build a truth set by requiring presence in some samples and absence
# in others.
# hmq and lmq can be either gatk or samtools formatted tables but must
# have exactly the same column layout.
# for samtools, ignore columns XXX
get.nonbulk.sites <- function(hmq, lmq, bulkgt, bulkref, bulkalt,
    bulk.min.dp=11, table.type='gatk')
{
    if (table.type == 'samtools')
        stop("samtools tables not supported at the moment")

    # strip data set of any site with evidence in bulk
    if (!missing(lmq)) {
        cat("using low mapping quality data. variant is only unsupported if neither low nor high mq reads allow it\n")
        # make lmq rows mirror hmq rows. lmq is NOT a proper superset of hmq
        lmq.aligned <- merge(hmq[,c('chr', 'pos', 'refnt', 'altnt')], lmq, all.x=T)

        newtab <- NULL
        if (table.type == 'gatk') {
            # gatk column format is repeated triples of (genotype, #ref, #alt)
            # first 7 columns are not sample data
            newtab <- hmq
            non.readcount.cols <- c(1:7,seq(8, ncol(hmq), 3))
            newtab[,-non.readcount.cols] <-
                pmax(hmq[,-non.readcount.cols], lmq.aligned[,-non.readcount.cols], na.rm=T)

            # sites missing in the lmq table might be poorly mappable and perhaps
            # should not be trusted
            newtab$missing.lmq <- is.na(lmq.aligned[,9])
        }
    } else {
        cat("no low mapping quality dataset supplied\n")
        newtab <- hmq
    }

    newtab <- newtab[newtab[,bulkref] + newtab[,bulkalt] >= bulk.min.dp &
                     newtab[,bulkalt] == 0 &
                     !grepl(",", newtab$altnt) &    # biallelic
                     newtab$dbsnp == '.',]

    if (table.type == 'gatk' & !missing(bulkgt))
        newtab <- newtab[newtab[,bulkgt] == '0/0',]

    cat("classifying sites..\n")
    k <- adapt.gatk(newtab)
    newtab$type <- k$type
    newtab$cells <- k$cells
    newtab$ncells <- k$ncells
    newtab$nout <- k$nout

    return(newtab)
}

# for 96 dimensional mut sigs
mutsig.cols <- rep(c('deepskyblue', 'black', 'firebrick2', 'grey', 'chartreuse3', 'pink2'), each=16)

get.3mer <- function(df) {
    require(BSgenome)
    require(BSgenome.Hsapiens.UCSC.hg19)

    comp <- c('A', 'C', 'G', 'T')
    names(comp) <- c('T', 'G', 'C', 'A')

    x <- df

    x$ctx <- getSeq(BSgenome.Hsapiens.UCSC.hg19,
                    names=paste("chr", x$chr, sep=''),
                    start=x$pos-1, end=x$pos+1, as.character=TRUE)
    x$ctx.rc <- sapply(strsplit(x$ctx, ""),
                    function(s) paste0(comp[s[c(3,2,1)]], collapse=''))

    x$type.and.ctx <- ifelse(x$refnt == 'C' | x$refnt == 'T',
                       paste0(x$ctx, ":", x$muttype),
                       paste0(x$ctx.rc, ":", x$muttype))
    x
}

plot.3mer <- function(x, ...) {
    bases <- c("A", 'C', 'G', 'T')

    # need to make a table of all possible contexts because they may not
    # be observed after filtering.
    t <- rep(0, 96)
    names(t) <- paste0(rep(bases, each=4),
                      rep(c('C', 'T'), each=48),
                      rep(bases, times=4),
                      ":",
                      rep(c('C', 'T'), each=48),
                      ">",
                      c(rep(c('A', 'G', 'T'), each=16),
                        rep(c('A', 'C', 'G'), each=16)))
    print(table(x$ctx))
    t2 <- table(x$type.and.ctx)
    t[names(t2)] <- t2
    tn <- do.call(rbind, strsplit(names(t), ":"))
    t <- t[order(tn[,2])]
    print(t)
    p <- barplot(t, las=3, col=mutsig.cols, names.arg=tn[order(tn[,2]), 1], space=0.5, border=NA, ...)
    abline(v=(p[seq(4,length(p)-1,4)] + p[seq(5,length(p),4)])/2, col='grey')
    legend('topright', ncol=2, legend=sort(unique(tn[,2])),
        fill=mutsig.cols[seq(1, length(mutsig.cols), 16)])
}

# to genotype doublets, run genotype.somatic twice: once in normal mode
# and again in doublet mode. then merge the results by unique sites.
# there will be differences in things like mda.pv, lysis.pv, as
# expected.
merge.doublet <- function(singlet, doublet) {
    x <- merge(singlet, doublet, by=colnames(doublet)[(1:101)[-53]],
          all=T, suffixes=c('', '.doublet'))
    x$subcell <- ifelse(is.na(x$lysis.pv), 'doublet specific',
        ifelse(is.na(x$lysis.pv.doublet), 'singlet specific', 'both'))
    x
}

# gatk is a dataframe containing GATK derived ref and alt counts
# genome-wide. it must contain the output of at least one single
# cell jointly called with at least one bulk.
# gatk.lowmq should contain the same data as 'gatk', but using a
#       very permissive mapping quality filter. This is to correct
#       for different read lengths between bulk and single cells.
#       IMPORTANT: the column layout of gatk and gatk.lowmq must match!
# sc.idx - column index of the single cell genotypes.  assumed
#          that columns sc.idx+1 and +2 are the ref and alt read
#          counts, respectively.
# bulk.idx - same as sc.idx, but for the jointly called bulk
# balprof.path - path to balance estimates. must contain files
#                named balprof.chrN.rda for N=1..22.
# chrs - only genotype chromosomes in this array. must be numeric
#        in {1..22}
#      NO LONGER WORKS
# call parameters
# bins - number of bins to use when determining artifact prevalence
# target.fdr - objective false discovery rate. no guarantee to
#              be matched.
# {sc,bulk}.min.dp - minimum read depth for calls
# sc.min.alt - minimum variant read support for calling sSNV
# bulk.max.alt - maximum variant read support tolerated in bulk
# max.gp.sd - sqrt(2) is where the transformation from a unimodal
#   distribution in normal space (where the gaussian process lives)
#   produces a bimodal distribution in allele fraction space [0,1].
#   for this reason, it is recommended to remove any larger SD.
#   removing the cutoff can improve sensitivity, but must
#   be used at your own risk.
#   be aware, however, that capping gp.sd can cause a large fraction
#   of the genome to be excluded (10-20%).
genotype.somatic <- function(gatk, gatk.lowmq, sc.idx, bulk.idx, somatic.ab,
    somatic.cigars, hsnp.cigars, cap.alpha=TRUE, doublet=FALSE, sites.only=FALSE,
    cg.id.q=0.9, cg.hs.q=0.9,
    random.seed=0, target.fdr=0.1, bins=20, max.gp.sd=Inf, #sqrt(2),
    bulkref=bulk.idx+1, bulkalt=bulk.idx+2, scref=sc.idx+1, scalt=sc.idx+2,
    bulk.max.alt=0, bulk.min.dp=11, sc.min.alt=2, sc.min.dp=6)
{
    call.fingerprint <- as.list(environment())
    # almost all of this information is saved in the results
    dont.save <- c('gatk', 'gatk.lowmq', 'somatic.ab', 'somatic.cigars')
    call.fingerprint <- call.fingerprint[!(names(call.fingerprint) %in% dont.save)]
    set.seed(random.seed)
    fraction.gp.sd.reliable <- NA
    if (sites.only) {
        cat("Skipping step 1, only returning somatic candidate and bulk sites.\n")
    } else {
        cat("step 1: finding het SNPs and candidate sSNVS..\n")
        gatk$muttype <- muttype.map[paste(gatk$refnt, gatk$altnt, sep=">")]
        gatk$dp <- gatk[,scalt] + gatk[,scref]
        gatk$af <- gatk[,scalt] / gatk$dp
        if (missing(somatic.cigars) | missing(hsnp.cigars))
            cat("    WARNING: no CIGAR data supplied. will increase FP rate\n")
    }

    # annotate each variant with the counts of supporting reads allowing
    # for low mapping quality.
    cat(sprintf("        attaching low MQ data..\n"))
    old.bulkref <- bulkref
    old.bulkalt <- bulkalt
    bulkref <- ncol(gatk) + bulkref - 2
    bulkalt <- ncol(gatk) + bulkalt - 2
    gatk <- merge(gatk, gatk.lowmq, by=c('chr', 'pos'),
        suffixes=c('', '.lmq'), all.x=TRUE)
    cat(sprintf("        adjusted bulk ref and alt columns to\n"))
    cat(sprintf("          -> ref=%s, alt=%s\n", colnames(gatk)[bulkref],
        colnames(gatk)[bulkalt]))
    cat(sprintf("        replacing %d NA low MQ counts with high MQ counts\n",
        sum(is.na(gatk[,bulkalt]))))
    gatk[,bulkref] <- ifelse(is.na(gatk[,bulkref]),
                             gatk[,old.bulkref], gatk[,bulkref])
    gatk[,bulkalt] <- ifelse(is.na(gatk[,bulkalt]),
                             gatk[,old.bulkalt], gatk[,bulkalt])


    cat("step 2: separating germline sites from somatic SNV candidates\n")
    cat(sprintf("        bulk min depth             %d\n", bulk.min.dp))
    cat(sprintf("        bulk max alt reads         %d\n", bulk.max.alt))
    cat(sprintf("        single cell min depth      %d\n", sc.min.dp))
    cat(sprintf("        single cell min alt reads  %d\n", sc.min.alt))
    germline <- gatk[gatk[,bulk.idx] == '0/1' &
                     gatk[,bulkref]+gatk[,bulkalt] >= bulk.min.dp,]
    somatic <- gatk[gatk[,bulkalt] <= bulk.max.alt &
                    gatk[,old.bulkalt] <= bulk.max.alt &
                    gatk[,bulkalt] + gatk[,bulkref] >= bulk.min.dp &
                    gatk[,scalt] >= sc.min.alt &
                    gatk[,scalt] + gatk[,scref] >= sc.min.dp &
                    # some sites with #alt=0 are called as variants presumably
                    # due to poor PairHMM likelihood ratios
                    gatk[,bulk.idx] == '0/0' &
                    # only retain biallelic sites
                    !grepl(",", gatk$altnt) &
                    # not present in dbSNP
                    gatk$dbsnp == '.',]

    cat(sprintf("        %d germline het SNVs\n", nrow(germline)))
    cat(sprintf("        %d somatic SNV candidates\n", nrow(somatic)))
    if (sites.only)
        return(list(germline=germline, somatic=somatic))

    cat(sprintf("        attaching AB estimates to somatic candidates..\n"))
    somatic <- merge(somatic, somatic.ab, all.x=T, by=c('chr', 'pos'))
    # choose the nearest AB
    if (doublet) {
        # XXX: might be true that i could get rid of the doublet factors
        # everywhere else in the analysis by setting AB=(closest AB)/2 here
        # not worrying about it now.
        ab1 <- 1/(1+exp(abs(somatic$gp.mu)))
        ab2 <- 1/(1+exp(-abs(somatic$gp.mu)))
        somatic$gp.mu <- ifelse(abs(somatic$af - ab1/2) < abs(somatic$af - ab2/2), abs(somatic$gp.mu), -abs(somatic$gp.mu))
        # choose the AB closest to the expected doublet AF
        somatic$ab <- ifelse(abs(somatic$af - ab1/2) < abs(somatic$af - ab2/2), ab1, ab2)
    } else {
        somatic$gp.mu <- ifelse(somatic$af < 1/2,
            -abs(somatic$gp.mu), abs(somatic$gp.mu))
        somatic$ab <- 1/(1+exp(-somatic$gp.mu))
    }
    if (max.gp.sd < Inf)
        cat("        removing sites with uncertain AB estimates\n")

    fraction.gp.sd.reliable <- mean(somatic$gp.sd < max.gp.sd)
    cat(sprintf("          -> %d / %d (%0.2f%%) sites retained\n",
        sum(somatic$gp.sd < max.gp.sd), nrow(somatic), 100*fraction.gp.sd.reliable))
    somatic <- somatic[somatic$gp.sd < max.gp.sd,]

    cat("step 3: computing p-values for filters\n")
    cat("        step 3a: allele balance consistency (somatic)\n")
    factor <- ifelse(doublet, 2, 1)
    somatic$abc.pv <-
        mapply(abc2, altreads=somatic[,scalt], gp.mu=somatic$gp.mu, gp.sd=somatic$gp.sd, factor=factor, dp=somatic$dp)

    cat("        step 3a: SKIPPING allele balance consistency (germline control)\n")
    # NB. includes NAs because of germline SNVs where SC depth=0
    #germline$abc.pv <-
        #mapply(abc, altreads=germline[,scalt], ab=germline$ab, dp=germline$dp)

    cat("        step 3b: lysis artifacts\n")
    factor <- ifelse(doublet, 4, 2)
    somatic$lysis.pv <-
        mapply(test2, altreads=somatic[,scalt], gp.mu=somatic$gp.mu, gp.sd=somatic$gp.sd, dp=somatic$dp, div=factor)

    cat("        step 3c: MDA artifacts\n")
    factor <- ifelse(doublet, 8, 4)
    somatic$mda.pv <-
        mapply(test2, altreads=somatic[,scalt], gp.mu=somatic$gp.mu, gp.sd=somatic$gp.sd, dp=somatic$dp, div=factor)


    cat(sprintf("step 4: controlling FDR = %0.3f\n", target.fdr))
    
    if (cap.alpha)
        cat(sprintf("        cap.alpha=TRUE: alpha <= %0.3f enforced despite artifact prevalence\n", target.fdr))

    somatic$popbin <- ceiling(somatic$af*bins)

    cat(sprintf("        estimating bounds on somatic mutation rate..\n"))
    max.dp <- quantile(germline$dp, prob=0.9)
    fcs <- lapply(sc.min.dp:max.dp, function(dp)
        fcontrol(germ.df=germline[germline$dp == dp & germline[,scalt] > sc.min.alt,],
                som.df=somatic[somatic$dp == dp,],
                bins=bins, doublet=doublet, min.dp=0)
    )
    fc.max <- fcontrol(germ.df=germline[germline$dp > max.dp & germline[,scalt] > 1,],
                som.df=somatic[somatic$dp > max.dp,],
                bins=bins, doublet=doublet, min.dp=0)
    fcs <- c(fcs, list(fc.max))

    burden <- as.integer(
        c(sum(sapply(fcs, function(fc) fc$est.somatic.burden[1])),  # min est
          sum(sapply(fcs, function(fc) fc$est.somatic.burden[2])))  # max est
    )

    cat(sprintf("        estimated callable somatic mutation burden range (%d, %d)\n",
        burden[1], burden[2]))
    cat("          -> using MAXIMUM burden\n")
    cat(sprintf("           --> WARNING: gp.sd not calculated for germline sites, no longer representative\n"))
    cat(sprintf("           --> warning: not estimating effect of ID filter\n"))
    

    cat("        step 4a: lysis artifact FDR\n")
    factor <- ifelse(doublet, 4, 2)
    nt.na <- mapply(function(dp, popbin)
        fcs[[min(dp, max.dp+1) - sc.min.dp + 1]]$pops$max[popbin,],
        somatic$dp, somatic$popbin)
    somatic <- cbind(somatic,
                lysis=t(mapply(match.fdr2,
                            gp.mu=somatic$gp.mu, gp.sd=somatic$gp.sd,
                            dp=somatic$dp,
                            nt=nt.na[1,],
                            na=nt.na[2,],
                            target.fdr=target.fdr, div=factor,
                            cap.alpha=cap.alpha)))

    cat("        step 4b: MDA artifact FDR\n")
    factor <- ifelse(doublet, 8, 4)
    somatic <- cbind(somatic,
                mda=t(mapply(match.fdr2,
                            gp.mu=somatic$gp.mu, gp.sd=somatic$gp.sd,
                            dp=somatic$dp,
                            nt=nt.na[1,],
                            na=nt.na[2,],
                            target.fdr=target.fdr, div=factor,
                            cap.alpha=cap.alpha)))

    cat("step 5: applying alignment filters (CIGAR.ID, CIGAR.HS)\n")
    if (missing(somatic.cigars) | missing(hsnp.cigars)) {
        cat("    WARNING: skipping alignment filtration. will increase FP rate\n")
    } else {
        somatic <- merge(somatic, somatic.cigars, by=c('chr', 'pos'), all.x=T)
        hsnp.cigars <- merge(hsnp.cigars, germline[,c('chr', 'pos', 'dp')], all.x=T)
        hsnp.cigars <- hsnp.cigars[hsnp.cigars$dp > 0,]
        somatic$cg.id.test <-
            test.against.germline(hsnp.cigars$ID/hsnp.cigars$dp,
                somatic=somatic$ID/somatic$dp, max.q=cg.id.q)
        somatic$cg.hs.test <-
            test.against.germline(hsnp.cigars$HS/hsnp.cigars$dp,
                somatic=somatic$HS/somatic$dp, max.q=cg.hs.q)
    }
    

    cat("step 6: calling somatic SNVs\n")
    somatic$pass <-
        somatic$abc.pv > 0.05 &
        somatic$lysis.pv <= somatic$lysis.alpha &
        somatic$mda.pv <= somatic$mda.alpha &
        somatic$cg.id.test & somatic$cg.hs.test
    cat(sprintf("        %d passing somatic SNVs\n", sum(somatic$pass)))
    cat(sprintf("        %d filtered somatic SNVs\n", sum(!somatic$pass)))

    # return non-data parameters and the output
    # the final RDA files are already large enough to be painful to work with
    return(c(call.fingerprint, burden=burden,
        list(fcs=fcs, somatic=somatic, germline=germline)))
}


# probability distribution of seeing y variant reads at a site with
# depth d with the estimate (gp.mu, gp.sd) of the local AB. model is
#     Y|p ~ Bin(d, 1/(1+exp(-b))),    b ~ N(gp.mu, gp.sd)
# the marginal distribution of Y (# variant reads) requires integrating
# over b, which has no closed form solution. we approximate it using
# gauss-hermite quadrature. increasing the number of nodes in ghd
# increases accuracy.
# 'factor' allows changing the relationship of ab -> af
# NOTE: for many calls, is best for the caller to compute ghd once
# and supply it to successive calls.
dreads <- function(ys, d, gp.mu, gp.sd, factor=1, ghd=gaussHermiteData(32)) {
    require(fastGHQuad)
    
    sapply(ys, function(y)
        ghQuad(function(x) {
                # ghQuad requires the weight function on x to be exp(-x^2)
                # (which is NOT a standard normal)
                b <- sqrt(2)*gp.sd*x + gp.mu
                exp(dbinom(y, size=d, prob=1/(factor*(1+exp(-b))), log=TRUE) - log(pi)/2)
            }, ghd
        )
    )
}


# determine the beta (power) as a function of alpha (FP rate)
# note that this is subject to heavy integer effects, specifically
# when D is low or when af is close to 0 or 1
# div controls the error model
estimate.alphabeta2 <- function(gp.mu, gp.sd, dp=30, alphas=10^seq(-4,0,0.05), div=2) {
    td <- data.frame(dp=0:dp,
        mut=dreads(0:dp, d=dp, gp.mu=gp.mu, gp.sd=gp.sd),
        err1=dreads(0:dp, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=div),
        err2=dreads(0:dp, d=dp, gp.mu=-gp.mu, gp.sd=gp.sd, factor=div))

    # (td[,4]+td[,5])/2 corresponds to an even mixture of errors on the p
    # and 1-p alleles.
    # XXX: i really don't know how I feel about this mixture. but use it
    # for now, maybe later come up with a better prior for the two error
    # modes.
    td$err <- (td$err1 + td$err2)/2
    td <- td[order(td$err),]
    td$cumerr <- cumsum(td$err)

    # FP probability assuming mutation occurs on this allele (wrt mu)
    fps <- sapply(alphas, function(alpha) {
        reject <- td$cumerr <= alpha
        # because there are two potential H0 models, a single alpha
        # does not make sense. Instead, use the maximum alpha to stay
        # conservative.
        max.alpha <- max(c(0, td$cumerr[td$cumerr <= alpha]))
        c(alpha=max.alpha, beta=sum(c(0, td$mut[reject == TRUE])))
    })

    list(td=td, alphas=fps[1,], betas=fps[2,], input.alphas=alphas)
}


# NOT equivalent to test(div=1)
abc2 <- function(altreads, gp.mu, gp.sd, dp, factor=1) {
    pmf <- dreads(0:dp, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=factor)
    p.value <- sum(pmf[pmf <= pmf[altreads + 1]])
    c(p.value=p.value)  # just names the return
}

test2 <- function(altreads, gp.mu, gp.sd, dp, div) {
    # equal mixture of errors in cis and trans
    err <- (dreads(0:dp, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=div) +
            dreads(0:dp, d=dp, gp.mu=-gp.mu, gp.sd=gp.sd, factor=div))/2

    p.value <- sum(err[err <= err[altreads + 1]])
    c(p.value=p.value)  # just names the return
}


# jitter - to remove integer effects. Using depth cutoffs is generally
# a bad strategy because DP and balance are correlated
bin.afs <- function(afs, bins=20, jit.amt=0.001) {
    sq <- seq(0, 1, 1/bins)
    jafs <- pmin(pmax(jitter(afs, amount=jit.amt), 0), 1)
    x <- findInterval(jafs, sq, left.open=TRUE)
    # af=0 will be assigned to bin 0 because intervals are (.,.]
    x[x==0] <- 1
    tx <- tabulate(x, nbins=bins)
    names(tx) <- apply(cbind(sq[-length(sq)], sq[-1]), 1, mean)
    tx
}


# the "rough" interval is not a confidence interval. it is just a
# heuristic that demonstrates the range of  reasonably consistent
# somatic mutation burdens
# WARNING! this is the *callable somatic burden*. if you want an
# estimate of the total somatic burden genome-wide, adjust this
# number by the fraction of genome represented in the somatic
# candidate set.
estimate.somatic.burden <- function(fc, min.s=1, max.s=5000, n.subpops=10, display=FALSE, rough.interval=0.99) {
    sim <- function(n.muts, g, s, n.samples=1000, diagnose=FALSE) {
        samples <- rmultinom(n=n.samples, size=n.muts, prob=g)
        if (diagnose) {
            boxplot(t(samples))
            lines(s, lwd=2, col=2)
        }
        mean(apply(samples, 2, function(col) all(col <= s)))
    }

    # determine how often a sample of N somatic mutations from the
    # proper het distribution (germlines) "fits" inside the somatic
    # candidate distribution.
    srange <- seq(min.s, max.s, length.out=n.subpops)
    fraction.embedded <- sapply(srange, sim, g=fc$g, s=fc$s)
    min.burden <- max(c(1, srange[fraction.embedded >= 1 - (1 - rough.interval)/2]))
    max.burden <- max(srange[fraction.embedded >= (1 - rough.interval)/2])
    c(min=min.burden, max=max.burden)
}


# {germ,som}.df need only have columns named dp and af
# estimate the population component of FDR
# ignore.100 - ignore variants with VAF=100
fcontrol <- function(germ.df, som.df, min.dp=10, jit.amt=0.0000001, bins=20, doublet=FALSE, rough.interval=0.99) {
    germ.afs <- germ.df$af[germ.df$dp >= min.dp & !is.na(germ.df$af)]
    som.afs <- som.df$af[som.df$dp >= min.dp & !is.na(som.df$af)]
    g <- bin.afs(germ.afs, jit.amt=jit.amt, bins=bins)  # counts, not probabilities
    s <- bin.afs(som.afs, jit.amt=jit.amt, bins=bins)

    # when fcontrol is used on small candidate sets (i.e., when controlling
    # for depth), there may be several 0 bins in s.
    g <- g*1*(s > 0)
    # XXX: i don't like this. it doesn't quite match the principle, but
    # rather addresses a limitation of the heuristic method.

    # returns (lower, upper) bound estimates
    approx.ns <- estimate.somatic.burden(fc=list(g=g, s=s),
        min.s=1, max.s=nrow(som.df), n.subpops=min(nrow(som.df), 100),
        rough.interval=rough.interval)

    pops <- lapply(approx.ns, function(n) {
        nt <- n*(g/sum(g))*1*(s > 0)
        # ensure na > 0, since FDR would be 0 for any alpha for na=0
        # XXX: the value 0.1 is totally arbitrary and might need to be
        # more carefully thought out.
        na <- pmax(s - nt, 0.1)
        cbind(nt=nt, na=na)
    })

    if (!doublet)
        return(list(est.somatic.burden=approx.ns,
             binmids=as.numeric(names(g)),
             g=g, s=s, pops=pops, doublet=doublet))

    # get the expected VAF distributions for UNSHARED mutations
    # (e.g., mean VAF=half of typical heterozygote
    # it is not correct to simply divide all AFs by 2.
    # use binomial sampling instead.
    gdf <- germ.df[germ.df$dp >= min.dp & !is.na(germ.df$af),]
    new.alt <- rbinom(n=nrow(gdf), size=round(gdf$af*gdf$dp), prob=1/2)
    germ.afs.unshared <- new.alt/gdf$dp
    g2 <- bin.afs(germ.afs.unshared, jit.amt=jit.amt, bins=bins)
    g2 <- g2*1*(s > 0)

    # remove the mutations already attributed to VAF=50% (i.e., shared)
    new.s <- pmax(s - pops$max[,1], 0)
    approx.ns2 <- estimate.somatic.burden(fc=list(g=g2, s=new.s),
        min.s=1, max.s=nrow(som.df), n.subpops=min(nrow(som.df), 100),
        rough.interval=rough.interval)

    pops <- lapply(approx.ns2, function(n) {
        nt <- n*(g2/sum(g2))*1*(new.s > 0)
        # ensure na > 0, since FDR would be 0 for any alpha for na=0
        na <- pmax(new.s - nt, 0.1)
        cbind(nt=nt, na=na)
    })

    list(est.somatic.burden=approx.ns2,
         binmids=as.numeric(names(g2)),
         g=g2, s=new.s, pops=pops, doublet=doublet)
}

# SUMMARY
# given a target FDR, try to match the FDR as best as we can without
# exceeding it.
#
# DETAILS
# for a given (af, dp), find the relation (alpha, beta) to determine
# what sensitivities and specificities are achievable.  using the
# estimated numbers of true mutations nt and artifacts na in the
# population, calculate the corresponding FDRs (alpha,beta,FDR).
# return (alpha, beta, FDR) such that
# FDR = arg max { FDR : FDR <= target.fdr }
# to break ties (when more than one (alpha,beta) produce the same
# FDR), select the triplet with minimal alpha and maximal beta.
# cap.alpha - because there is often such a strong spike at low VAF
#             amongst somatic candidate sites, the FDR control
#             procedure sometimes results in prior error rates of <1%
#             amongst high VAF candidate sSNVs. In this case, it will
#             often accept variants that are very clearly more consistent
#             with the artifact model(s) than with a true mutation. The
#             cap.alpha parameter helps to address this issue by requiring
#             alpha <= target FDR.
#             To be more clear: when nt >> na, FDR matching can pass
#             variants that have very high alpha (e.g., cases with
#             alpha ~ 0.4, beta ~ 0.7). This puts significant pressure
#             on the accuracy of the nt and na estimates.
match.fdr2 <- function(gp.mu, gp.sd, dp, nt, na, target.fdr=0.1, div=2, cap.alpha=TRUE) {
    alphabeta <- estimate.alphabeta2(gp.mu=gp.mu, gp.sd=gp.sd, dp=dp, div=div)
    x <- data.frame(alpha=alphabeta$alphas, beta=alphabeta$betas,
               fdr=ifelse(alphabeta$alphas*na + alphabeta$betas*nt > 0,
                          alphabeta$alphas*na / (alphabeta$alphas*na + alphabeta$betas*nt), 0))
    x <- rbind(c(0, 0, 0), x)  # when no parameter meet target fdr
    # use (a,) from the highest FDR <= target.fdr
    x <- x[x$fdr <= target.fdr,] 
    if (cap.alpha)
        x <- x[x$alpha <= target.fdr,]
    x <- x[x$fdr == max(x$fdr),]  # may be more than one row
    # breaking ties by minimum alpha, maximum beta
    x <- x[x$alpha == min(x$alpha),]
    x <- x[x$beta == max(x$beta),]
    unlist(x[1,])
}




###############################################################################
# various plotting routines
###############################################################################


plot.ab <- function(ab) {
    layout(matrix(1:4,nrow=2,byrow=T))
    td <- ab$td[order(ab$td$dp),]
    plot(td$dp, td$mut, type='l', ylim=range(td[,c('mut', 'err1', 'err2')]),
        xlab="Depth", ylab="Model probabiblity")
    lines(td$dp, pmax(td$err1,td$err2), lty='dotted', col=2)
    plot(x=ab$input.alphas, y=ab$alphas, xlab="Requested alpha", ylab="Estimated alpha", log="xy")
    plot(ab$alphas, ab$betas, xlab='FP rate', ylab='Power')
    abline(h=0, lty='dotted')
    plot(ab$alphas, ab$betas, log='x', xlab='log(FP rate)', ylab='Power')
    abline(h=0, lty='dotted')
}


plot.3mer <- function(x, ...) {
    bases <- c("A", 'C', 'G', 'T')

    # need to make a table of all possible contexts because they may not
    # be observed after filtering.
    t <- rep(0, 96)
    names(t) <- paste0(rep(bases, each=4),
                      rep(c('C', 'T'), each=48),
                      rep(bases, times=4),
                      ":",
                      rep(c('C', 'T'), each=48),
                      ">",
                      c(rep(c('A', 'G', 'T'), each=16),
                        rep(c('A', 'C', 'G'), each=16)))
    print(table(x$ctx))
    t2 <- table(x$type.and.ctx)
    t[names(t2)] <- t2
    tn <- do.call(rbind, strsplit(names(t), ":"))
    t <- t[order(tn[,2])]
    print(t)
    p <- barplot(t, las=3, col=mutsig.cols, names.arg=tn[order(tn[,2]), 1], space=0.5, border=NA, ...)
    abline(v=(p[seq(4,length(p)-1,4)] + p[seq(5,length(p),4)])/2, col='grey')
    legend('topright', ncol=2, legend=sort(unique(tn[,2])),
        fill=mutsig.cols[seq(1, length(mutsig.cols), 16)])
}


plot.fcontrol <- function(fc) {
    layout(matrix(1:(1+length(fc$pops)), nrow=1))
    plot(fc$binmids, fc$g/sum(fc$g),
        ylim=range(c(fc$g/sum(fc$g), fc$s/sum(fc$s))),
        type='l', lwd=2)
    lines(fc$binmids, fc$s/sum(fc$s), col=2, lwd=2)

    for (i in 1:length(fc$pops)) {
        pop <- fc$pops[[i]]
        barplot(names.arg=fc$binmids, t(pop), col=1:2,
            main=sprintf('Assumption: ~%d true sSNVs', sum(round(pop[,1],0))), las=2)
        legend('topright', fill=1:2, legend=c('Ntrue', 'Nartifact'))
    }
}


# using the (alpha, beta) relationships and fcontrol population
# estimations, determine average sensitivity per AF with a
# (theoretically) controlled FDR
plot.fdr <- function(fc, dps=c(10,20,30,60,100,200), target.fdr=0.1, div=2) {
    afs <- fc$binmids
    layout(matrix(1:(3*length(fc$pops)), nrow=3))
    for (i in 1:length(fc$pops)) {
        # from fcontrol: pops[[i]] has rows corresponding to afs
        pop <- fc$pops[[i]]
        l <- lapply(dps, function(dp) {
            sapply(1:length(afs), function(i)
                match.fdr(afs[i], dp, nt=pop[i,1], na=pop[i,2],
                    target.fdr=target.fdr, div=div)
            )
        })
        matplot(x=afs, sapply(l, function(ll) ll[3,]), type='l', lty=1,
            main=sprintf("Assuming %d true sSNVs", sum(pop[,1])),
            xlab="AF (binned)", ylab="FDR", ylim=c(0, 1.1*target.fdr))
        abline(h=target.fdr, lty='dotted')
        matplot(x=afs, sapply(l, function(ll) ll[1,]), type='l', lty=1,
            xlab="AF (binned)", ylab="log(alpha)", log='y', ylim=c(1e-5,1))
        abline(h=10^-(1:5), lty=2)
        matplot(x=afs, sapply(l, function(ll) ll[2,]), type='l', lty=1,
            xlab="AF (binned)", ylab="Power", ylim=0:1)
    }
}


# must have a column named 'cells', produced by cell.inclusion above
# or can be specified by 'colname'
plot.cell.sharing <- function(df, nsamples=13, colname='cells') {
    require(UpSetR)
    upset(fromExpression(table(df[[colname]])), nsets=nsamples, order.by='freq')
}


plot.ssnv.region <- function(chr, pos, alt, ref, fits, fit.data, upstream=5e4, downstream=5e4,
    # needed for infer.gp
    source("predict.R")

    gp.extend=1e5, n.gp.points=100) {
    d <- fit.data[fit.data$chr==chr & fit.data$pos >= pos - upstream &
                  fit.data$pos <= pos + downstream,]

    cat("estimating AB in region..\n")
    est.at <- seq(pos - upstream, pos + downstream, length.out=n.gp.points)
    fit.chr <- as.data.frame(fits[chr,,drop=FALSE])
    colnames(fit.chr) <- c('a', 'b', 'c', 'd', 'logq')
    gp <- infer.gp(vars=data.frame(pos=est.at),
        fit=fit.chr, chr=fit.data[fit.data$chr == chr,],
        chunk=250, flank=gp.extend)
    gp$pos <- est.at
    plot.gp.confidence(df=gp, add=FALSE)
    points(d$pos, d$hap1/d$dp, pch=20, ylim=0:1)
    points(pos, alt/(alt+ref), pch=20, cex=1.25, col=2)
}


# NOTE: the 95% probability interval is in the NORMAL space, not the
# fraction space [0,1]. After the logistic transform, the region in
# fraction space may not be an HDR.
plot.gp.confidence <- function(pos, gp.mu, gp.sd, df, sd.mult=2,
    logspace=FALSE, tube.col=rgb(0.9, 0.9, 0.9),
    line.col='black', tube.lty='solid', line.lty='solid', add=TRUE)
{
    if (!missing(df)) {
        pos <- df$pos
        gp.mu <- df$gp.mu
        gp.sd <- df$gp.sd
    }

    # maybe the analytical solution exists, but why derive it
    if (!logspace) {
        cat("transforming to AF space...\n")
        sd.upper <- 1 / (1 + exp(-(gp.mu + sd.mult*gp.sd)))
        sd.lower <- 1 / (1 + exp(-(gp.mu - sd.mult*gp.sd)))
        gp.mu <- 1 / (1 + exp(-gp.mu))
    } else {
        sd.upper <- gp.mu + sd.mult*gp.sd
        sd.lower <- gp.mu - sd.mult*gp.sd
    }

    if (!add) {
        plot(NA, NA, xlim=range(pos), ylim=0:1)
    }

    polygon(c(pos, rev(pos)), c(gp.mu, rev(sd.lower)),
        col=tube.col, border=line.col, lty=tube.lty)
    polygon(c(pos, rev(pos)), c(gp.mu, rev(sd.upper)),
        col=tube.col, border=line.col, lty=tube.lty)
    lines(pos, gp.mu, lwd=2, col=line.col, lty=line.lty)
}
