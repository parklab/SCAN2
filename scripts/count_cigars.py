#!/usr/bin/env python

import pysam
import argparse
import os.path

ap = argparse.ArgumentParser(
    "Count CIGAR operations for all reads overlapping a given position. "
    "This does not count CIGAR operations overlapping _the position_, "
    "which means that it will notice if, for example, there is a large "
    "cluster of clipping operations nearby (i.e., within roughly the "
    "size of a read).")
ap.add_argument('input_bam')
ap.add_argument('input_list')
ap.add_argument('output_txt')
args = ap.parse_args()

if os.path.isfile(args.output_txt):
    raise RuntimeError('output file "%s" already exists, please remove it first' % args.output_txt)

samfile = pysam.AlignmentFile(args.input_bam, "rb")
with open(args.output_txt, 'w') as outfile:
    outfile.write('chr\tpos\tM.cigars\tID.cigars\tHS.cigars\tother.cigars\tdp.cigars\n')
    with open(args.input_list, 'r') as posfile:
        for line in posfile:
            # header line
            if line.startswith('chr\tpos\t'):
                continue
    
            chrom, pos, refnt, altnt = line.strip().split("\t")
            pos = int(pos)

            match = 0
            indel = 0
            clip = 0
            other = 0
            dp = 0
            reads = samfile.fetch(chrom, pos, pos+1)
            for r in reads:
                if r.mapping_quality >= 60 and r.is_paired and \
                    not r.is_duplicate and not r.is_qcfail and \
                    not r.is_secondary and not r.is_supplementary:
                    # e.g., 30M1I -> [(CMATCH,30), (CINS,1), ...]
                    tups = r.cigartuples
                    cops = [ tup[0] for tup in tups ]
                    this_match = cops.count(pysam.CMATCH)
                    this_indel = cops.count(pysam.CINS) + cops.count(pysam.CDEL)
                    this_clip = cops.count(pysam.CHARD_CLIP) + cops.count(pysam.CSOFT_CLIP)
                    this_other = len(cops) - (this_match + this_indel + this_clip)
                    match = match + this_match
                    indel = indel + this_indel
                    clip = clip + this_clip
                    other = other + this_other
                    dp = dp + 1

            outfile.write(chrom + '\t' + str(pos) + '\t' + str(match) + '\t' +
                str(indel) + '\t' + str(clip) + '\t' + str(other) + '\t' +
                str(dp) + '\n')
