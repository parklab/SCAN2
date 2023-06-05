#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
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

# `data` contains 1kb-resolution bins tiling the whole genome. This table
# is often ~1GB in size.
results@spatial.sensitivity$data <- NULL

# R's glm saves a copy of all data used to fit a model, so we keep only the model
# summaries in the minimized object.  However, even the model summaries contain
# a copy of the entire environment used for fitting in an attribute named
#      attr(summary(model)$terms, '.Environment')
# It appears that deleting this environment does not affect displaying basic info
# about the model, such as fitted coefs and p-values.
results@spatial.sensitivity$models <- lapply(results@spatial.sensitivity$models, function(m) {
    s <- summary(m)
    attr(s$terms, '.Environment') <- c()
    s
})

save(results, file=out.rda)
