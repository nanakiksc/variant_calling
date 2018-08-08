#!/usr/bin/env bash

SM=$1

HOME=/home/pcusco
BWANC=4
VQUAL=30
BUILD=hg19
DATAD=${HOME}/data/${BUILD}
FASTA=${DATAD}/masked/${BUILD}.karyosorted.fa
BWAIX=${DATAD}/masked/bwa/${BUILD}.masked
DBSNP=${DATAD}/dbsnp/All_20180423_chr.vcf.gz # chrsorted
MILLS=${DATAD}/gatk_resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
PHASE=${DATAD}/gatk_resources/1000G_phase1.indels.hg19.sites.vcf.gz
CPTUR=${HOME}/avivancos/plasma_umis/capture_regions/171106_HG19_VHIO_UMIs_EZ_HX3_capture_targets.bed
ANNDB=${HOME}/utils/annovar/humandb/

#bwa mem -t ${BWANC} -M -R '@RG\tID:1\tLB:1\tPL:illumina\tSM:1\tPU:1.1.1' \
#    ${BWAIX} \
#    <(zcat ${SM}_*_R1_*.fastq.gz) \
#    <(zcat ${SM}_*_R2_*.fastq.gz) | \
#    samtools sort -o ${SM}.raw.bam

#if [ $? -ne 0 ]; then exit $?; fi

#gatk MarkDuplicatesWithMateCigar \
#    --INPUT ${SM}.raw.bam \
#    --METRICS_FILE ${SM}.metrics.txt \
#    --OUTPUT ${SM}.dedup.bam

#if [ $? -ne 0 ]; then exit $?; fi
#rm ${SM}.raw.bam

gatk ReorderSam \
    --INPUT ${SM}.dedup.bam \
    --OUTPUT ${SM}.reorder.bam \
    --REFERENCE ${FASTA}

if [ $? -ne 0 ]; then exit $?; fi
#rm ${SM}.dedup.bam

samtools index ${SM}.reorder.bam

if [ $? -ne 0  ]; then exit $?; fi

java -jar ${HOME}/utils/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
    --analysis_type RealignerTargetCreator \
    --reference_sequence ${FASTA} \
    --input_file ${SM}.reorder.bam \
    --known ${MILLS} \
    --known ${PHASE} \
    --out ${SM}.intervals \
    --filter_bases_not_stored \
    --intervals ${CPTUR}

if [ $? -ne 0 ]; then exit $?; fi

java -jar ${HOME}/utils/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
    --analysis_type IndelRealigner \
    --reference_sequence ${FASTA} \
    --input_file ${SM}.reorder.bam \
    --targetIntervals ${SM}.intervals \
    --out ${SM}.realign.bam \
    --knownAlleles ${MILLS} \
    --knownAlleles ${PHASE} \
    --maxReadsForRealignment 150000 \
    --intervals ${CPTUR} \
    --filter_bases_not_stored

if [ $? -ne 0 ]; then exit $?; fi
#rm ${SM}.reorder.bam ${SM}.reorder.bai ${SM}.intervals

gatk BaseRecalibrator \
    --input ${SM}.realign.bam \
    --known-sites ${DBSNP} \
    --known-sites ${MILLS} \
    --known-sites ${PHASE} \
    --output ${SM}.recal.table \
    --reference ${FASTA}

if [ $? -ne 0 ]; then exit $?; fi

gatk ApplyBQSR \
    --bqsr-recal-file ${SM}.recal.table \
    --input ${SM}.realign.bam \
    --output ${SM}.recal.bam

if [ $? -ne 0 ]; then exit $?; fi
#rm ${SM}.realign.bam ${SM}.recal.table

gatk HaplotypeCaller \
    --input ${SM}.recal.bam \
    --output ${SM}.vcf \
    --reference ${FASTA} \
    --dbsnp ${DBSNP} \
    --standard-min-confidence-threshold-for-calling ${VQUAL} \
    --intervals ${CPTUR}

if [ $? -ne 0 ]; then exit $?; fi

${HOME}/utils/annovar/convert2annovar.pl \
    --format vcf4 \
    --includeinfo \
    --outfile ${SM}.ann \
    ${SM}.vcf

if [ $? -ne 0 ]; then exit $?; fi
#rm ${SM}.vcf

${HOME}/utils/annovar/annotate_variation.pl \
    --outfile ${SM} \
    --buildver ${BUILD} \
    --hgvs \
    ${SM}.ann \
    ${ANNDB}

if [ $? -ne 0 ]; then exit $?; fi
#rm ${SM}.ann

${HOME}/src/Pipelines/variant_calling/format.py ${SM}.exonic_variant_function > ${SM}.xls

if [ $? -ne 0 ]; then exit $?; fi
#rm ${SM}.*variant_function ${SM}.log
exit 0
