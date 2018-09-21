#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys

print '\t'.join(['dbSNP', 'Locus', 'Type', 'Transcripts', 'Ref', 'Alt', 'Reads1', 'Reads2', 'VarFreq'])

exonic_variant_function = sys.argv[1] + '.exonic_variant_function'
variant_function = sys.argv[1] + '.variant_function'

def parse_variants(line, offset):
    var_type = sline[0 + offset]
    dbsnp = sline[9 + offset] if sline[9 + offset] != '.' else ''
    chrom = sline[2 + offset].lstrip('chr')
    start = sline[3 + offset]
    end = sline[4 + offset]
    locus = '%s:%s-%s' % (chrom, start, end)
    transcript = sline[1 + offset].rstrip(',')
    ref = sline[5 + offset]
    alt = sline[6 + offset]
    format_field = sline[16 + offset].split(':')
    read1 = format_field[10]
    read2 = format_field[2]
    varfreq = format_field[6]

    return '\t'.join([dbsnp, locus, var_type, transcript, ref, alt, read1, read2, varfreq])


with open(exonic_variant_function) as fin:
    for line in fin:
        sline = line.split('\t')

        var_type = sline[1]
        if var_type == 'synonymous SNV':
            continue

        print parse_variants(sline, 1)

with open(variant_function) as fin:
    for line in fin:
        sline = line.split('\t')

        var_type = sline[0]
        if var_type == 'splicing':
            print parse_variants(sline, 0)

