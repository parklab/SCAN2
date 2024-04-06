#!/bin/bash

if [ $# -lt 3 ]; then
    echo 'usage: concat_tabix.sh "bed"|"vcf" output.tab.gz in1.tab.gz [ in2.tab.gz ... inN.tab.gz ]'
    echo "the first argument (bed|vcf) is the argument to the tabix preset (-p) index option"
    exit 1
fi

preset=$1

if [ "x$preset" != "xvcf" ] && [ "x$preset" != "xbed" ]; then
    echo "argument 1 must be either 'vcf' or 'bed', case sensitive"
    exit 1
fi

outfile=$2

shift 2     # now $@ is all input files
infile1=$1

if [ -f $outfile ]; then
    echo "output file $outfile already exists, please delete it and try again"
    exit 1
fi

echo "concatenating $# input files.."
(tabix --only-header $infile1 ;
    for infile in $@; do gunzip -c $infile | grep -v '^#'; done) \
| bgzip -c > $outfile

echo "indexing $outfile.."
tabix -p $preset $outfile
