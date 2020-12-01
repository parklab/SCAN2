sc <- read.table(snakemake@input[['sc']], header=T, stringsAsFactors=F)
bulk <- read.table(snakemake@input[['bulk']], header=T, stringsAsFactors=F)
cd <- merge(sc, bulk, by=c('chr', 'pos'), suffixes=c('', '.bulk'))

cigar.emp.score <- function(training, test, which=c('id', 'hs')) {
    xt <- training[,paste0(which, '.score.x')]
    yt <- training[,paste0(which, '.score.y')]
    x <- test[,paste0(which, '.score.x')]
    y <- test[,paste0(which, '.score.y')]
    mapply(function(xi, yi) mean(xt >= xi & yt >= yi, na.rm=T), x, y)
}

cd$id.score.y <- cd$ID.cigars / cd$dp.cigars
cd$id.score.x <- cd$ID.cigars.bulk / cd$dp.cigars.bulk
cd$id.score <- cigar.emp.score(training=cd, test=cd, which='id')
cd$hs.score.y <- cd$HS.cigars / cd$dp.cigars
cd$hs.score.x <- cd$HS.cigars.bulk / cd$dp.cigars.bulk
cd$hs.score <- cigar.emp.score(training=cd, test=cd, which='hs')

cigar.training <- cd
save(cigar.emp.score, cigar.training, file=snakemake@output[[1]])
