#!/usr/bin/env Rscript
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 4:00:00
#SBATCH --mem=24G

# Simple wrapper of the enrichment.r code to analyze enrichment for
# any genomic intervals defined in a BED file.
source("/n/data1/hms/dbmi/park/jluquette/glia/scripts/enrichment.r")

command.line.analysis(function(bedf)
    enrich.data(gbed=read.inout.bed(bedf, add.chr.prefix=FALSE)))
