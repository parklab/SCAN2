#!/usr/bin/env Rscript
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 2:00:00
#SBATCH --mem=32G

# Simple wrapper of the enrichment.r code to analyze enrichment for
# any genomic intervals defined in a BED file.
source("/n/data1/hms/dbmi/park/jluquette/glia/scripts/enrichment.r")

command.line.analysis(function(bedf)
    enrich.data(gbed=read.feature.bed(bedf, feature.name='bed.feature', remove.chr.prefix=TRUE)))

print(gc())
