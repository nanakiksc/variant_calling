#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys

print '\t'.join(['dbSNP', 'Locus', 'Type', 'Transcripts', 'Ref', 'Alt', 'Reads1', 'Reads2', 'VarFreq', 'Pvalue', 'PolyPhen2_score', 'PolyPhen2_pred', 'MAF_eur', 'ClinVar'])

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
    read1 = format_field[4]
    read2 = format_field[5]
    varfreq = format_field[6]
    pval = format_field[7]
    info_field = sline[14 + offset].split(';')
    info_dict = {key: value for key, value in (field.split('=') for field in info_field)}
    poly_score = info_dict['Polyphen2_HVAR_score'] if 'Polyphen2_HVAR_score' in info_dict else ''
    poly_pred = info_dict['Polyphen2_HVAR_pred'] if 'Polyphen2_HVAR_pred' in info_dict else ''
    maf = info_dict['AF_NFE'] if 'AF_NFE' in info_dict else ''
    clinvar = info_dict['CLNSIG'] if 'CLNSIG' in info_dict else ''

    return '\t'.join([dbsnp, locus, var_type, transcript, ref, alt, read1, read2, varfreq, pval, poly_score, poly_pred, maf, clinvar])


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

