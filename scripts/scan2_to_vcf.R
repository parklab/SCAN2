#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3)
    stop("usage: scan2_to_vcf.R in.rda out.vcf {simple.filters:true|false}")

in.rda <- args[1]
out.file <- args[2]
simple.filters <- as.logical(args[3])

if (file.exists(out.file) & out.file != '/dev/stdout')
    stop(sprintf("output file %s already exists, please delete it first", out.file))

if (out.file == '/dev/stdout') {
    devnull <- file('/dev/null', 'w')
    sink(file=devnull, type='message')
}

library(scan2)

message('loading SCAN2 object..')
object <- get(load(in.rda))
message("writing to ", out.file)
write.vcf(object=object, file=out.file, simple.filters=simple.filters)
    # hack for older objects without @config
    #config=list(ref='/n/data1/hms/dbmi/park/jluquette/glia/analysis/resources/human_g1k_v37_decoy.fasta'))
