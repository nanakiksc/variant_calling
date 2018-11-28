#!/usr/bin/env bash

SM=$1
MINAD=$2

HOME=/home/pcusco
BUILD=hg19
DATAD=${HOME}/data/${BUILD}
FASTA=${DATAD}/masked/${BUILD}.chr_sorted.fa
CPTUR=${HOME}/avivancos/plasma_umis/capture_regions/171106_HG19_VHIO_UMIs_EZ_HX3_capture_targets.bed

samtools mpileup -B \
    -f ${FASTA} \
    -o ${SM}.mpileup \
    ${SM}.recal.bam

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#rm ${SM}.recal.bam

samtools depth -a \
    -b ${CPTUR} \
    ${SM}.recal.bam | \
    awk '{ n += 1; s += $3 } END { printf "%.0f", s / n / 100 }' \
    > ${SM}.panel_avg_cov.txt

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

MINCV=$(cat ${SM}.panel_avg_cov.txt)
MINCV=$(( ${MINCV} > ${MINAD} ? ${MINCV} : ${MINAD} ))

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

java -jar ${HOME}/utils/varscan/VarScan-2.4.x/VarScan.v2.4.3.jar \
    mpileup2cns ${SM}.mpileup \
    --min-coverage ${MINCV} \
    --min-reads2 ${MINAD} \
    --min-var-freq 0.0 \
    --p-value 0.05 \
    --output-vcf 1 \
    --variants | \
    bgzip > ${SM}.vcf.gz

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.mpileup ${SM}.panel_avg_cov.txt

tabix ${SM}.vcf.gz

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
exit 0
