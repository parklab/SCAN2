#!/usr/bin/env python

import re
import argparse
import sys

def count_line(cigars):
    if cigars.strip() == '':
        return (0, 0, 0, 0)

    # reduce just to cigar operations
    s = ''.join([ re.sub('[0-9\ ]', '', cigar) for cigar in cigars ])
    m = s.count('M')
    indel = s.count('I') + s.count('D')
    clips = s.count('S') + s.count('H')
    other = len(s) - m - indel - clips
    return (m, indel, clips, other)


ap = argparse.ArgumentParser()
ap.add_argument("raw_cigars")
args = ap.parse_args()

#with open('raw_germline_cigars.txt', 'r') as f:
firstline = True
with open(args.raw_cigars, 'r') as f:
    for line in f:
        if firstline:
            samples = line.strip().split("\t")
            header = 'chr\tpos'
            for s in samples[2:]:
                header = header + '\t' + \
                    '\t'.join([ x + '.' + s for x in [ 'M', 'ID', 'HS', 'other' ] ])
            print(header)
            firstline = False
            continue
        line = line.strip()
        fields = line.split('\t')
        chrom = fields[0]
        pos = fields[1]


        sys.stdout.write(chrom + "\t" + pos)
        if len(fields) < 3:
            sys.stdout.write('\t' + '\t'.join([ str(x) for x in count_line('') ]))

        for cigars in fields[2:]:
            sys.stdout.write('\t' + '\t'.join([ str(x) for x in count_line(cigars) ]))
        sys.stdout.write('\n')
