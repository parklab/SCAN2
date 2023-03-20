#!/usr/bin/env Rscript

if (length(args) != 2)
    stop("usage: script.R in.rda out.rda")

in.rda <- args[1]
out.rda <- args[2]

if (file.exists(out.rda)) 
    stop(paste("output file", out.rda, "already exists, please delete it first"))


library(scan2)

load(in.rda)

if ('rescue' %in% colnames(results@gatk)) {
    results@gatk <- results@gatk[(pass == TRUE | rescue == TRUE) & somatic.candidate==TRUE]
} else {
    results@gatk <- results@gatk[pass == TRUE & somatic.candidate==TRUE]
}

save(results, file=out.rda)
