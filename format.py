#!/usr/bin/env python
#-*- coding:utf-8 -*-

# Hard filtering rules.
MIN_COV = 8
MIN_ALT = 7

import sys

print '\t'.join(['dbSNP', 'Locus', 'Type', 'Transcripts', 'Ref', 'Alt', 'Reads1', 'Reads2', 'VarFreq'])

with open(sys.argv[1]) as fin:
    for line in fin:
        sline = line.split('\t')

        var_type = sline[1]
        if var_type == 'synonymous SNV':
            continue
        dbsnp = sline[10] if sline[10] != '.' else ''
        chrom = sline[3].lstrip('chr')
        start = sline[4]
        end = sline[5]
        locus = '%s:%s-%s' % (chrom, start, end)
        transcript = sline[2].rstrip(',')
        ref = sline[6]
        alt = sline[7]
        format_field = sline[17].split(':')
        read1, read2 = format_field[1].split(',')
        depth = format_field[2]
        assert int(read1) + int(read2) == int(depth)
#        if int(depth) < MIN_COV or int(read2) < MIN_ALT:
#            continue
        varfreq = str(round(float(read2) / float(depth) * 100, 2)).replace('.', ',') + '%'

        print '\t'.join([dbsnp, locus, var_type, transcript, ref, alt, read1, read2, varfreq])
