#!/usr/bin/env python
#-*- coding:utf-8 -*-

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
        read1, read2 = format_field[4:6]
        varfreq = format_field[6]

        print '\t'.join([dbsnp, locus, var_type, transcript, ref, alt, read1, read2, varfreq])
