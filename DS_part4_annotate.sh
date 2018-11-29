#!/usr/bin/env bash

SM=$1

HOME=/home/pcusco
BUILD=hg19
DATAD=${HOME}/data/${BUILD}
DBSNP=${DATAD}/dbsnp/All_20180423_chr.vcf.gz
CPTUR=${HOME}/avivancos/plasma_umis/capture_regions/171106_HG19_VHIO_UMIs_EZ_HX3_capture_targets.bed
CLNVR=${DATAD}/ClinVar/clinvar_20180701_chr.vcf.gz
GNMAD=${DATAD}/gnomAD/exomes/gnomad.exomes.r2.0.2.sites_chr.vcf.gz
DNSFP=${DATAD}/dbnsfp33a/hg19_dbnsfp33a.vcf.gz
ANNDB=${HOME}/utils/annovar/humandb/

bcftools annotate \
    --annotations ${DBSNP} \
    --columns ID \
    --regions-file ${CPTUR} \
    ${SM}.vcf.gz | \
    bgzip > ${SM}.dbsnp.vcf.gz

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.vcf.gz*

tabix ${SM}.dbsnp.vcf.gz

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

bcftools annotate \
     --annotations ${CLNVR} \
     --columns INFO/CLNSIG \
     --regions-file ${CPTUR} \
     ${SM}.dbsnp.vcf.gz | \
     bgzip > ${SM}.clinvar.vcf.gz

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.dbsnp.vcf.gz*

tabix ${SM}.clinvar.vcf.gz

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

bcftools annotate \
     --annotations ${GNMAD} \
     --columns INFO/AF_NFE \
     --regions-file ${CPTUR} \
     ${SM}.clinvar.vcf.gz | \
     bgzip > ${SM}.gnomad.vcf.gz

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.clinvar.vcf.gz*

tabix ${SM}.gnomad.vcf.gz

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

bcftools annotate \
     --annotations ${DNSFP} \
     --columns INFO/Polyphen2_HVAR_score,INFO/Polyphen2_HVAR_pred \
     --output ${SM}.poly.vcf \
     --regions-file ${CPTUR} \
     ${SM}.gnomad.vcf.gz

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.gnomad.vcf.gz*

${HOME}/utils/annovar/convert2annovar.pl \
    --format vcf4 \
    --includeinfo \
    --outfile ${SM}.ann \
    ${SM}.poly.vcf

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.poly.vcf

${HOME}/utils/annovar/annotate_variation.pl \
    --outfile ${SM} \
    --buildver ${BUILD} \
    --hgvs \
    ${SM}.ann \
    ${ANNDB}

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.ann

./format.py ${SM} > ${SM}.xls

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.*variant_function ${SM}.log
exit 0
