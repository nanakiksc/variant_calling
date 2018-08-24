#!/usr/bin/env bash

#I1=$1
#I2=$2
SM=$1

HOME=/home/pcusco
BWANC=4
VQUAL=30
BUILD=hg19
DATAD=${HOME}/data/${BUILD}
FASTA=${DATAD}/masked/${BUILD}.karyosorted.fa
BWAIX=${DATAD}/masked/bwa/${BUILD}.masked
DBSNP=${DATAD}/dbsnp/All_20180423_chr.vcf.gz # chr
MILLS=${DATAD}/gatk_resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
PHASE=${DATAD}/gatk_resources/1000G_phase1.indels.hg19.sites.vcf.gz
CPTUR=${HOME}/avivancos/plasma_umis/capture_regions/171106_HG19_VHIO_UMIs_EZ_HX3_capture_targets.bed
ANNDB=${HOME}/utils/annovar/humandb/

#fastq_quality_filter -q 30 -p 75 \
#    -i <(zcat ${I1}) \
#    -o ${SM}_filtered_R1.fastq
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#
#fastq_quality_filter -q 30 -p 75 \
#    -i <(zcat ${I2}) \
#    -o ${SM}_filtered_R2.fastq
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#
#${HOME}/utils/Homer/bin/homerTools trim \
#    -5 TGACT \
#    -minMatchLength 5 \
#    -min 55 \
#    ${SM}_filtered_R1.fastq \
#    ${SM}_filtered_R2.fastq
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#rm ${SM}_filtered_R*.fastq *.lengths
#
#${HOME}/utils/Homer/bin/homerTools trim \
#    -3 AGATCGGAAGAGCACACGTCT \
#    -mis 2 \
#    -minMatchLength 4 \
#    -min 55 \
#    ${SM}_filtered_R1.fastq.trimmed \
#    ${SM}_filtered_R2.fastq.trimmed
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#rm ${SM}_filtered_R*.fastq.trimmed *.lengths
#
#${HOME}/utils/pairfq/Pairfq-0.17.0/bin/pairfq addinfo \
#    -i ${SM}_filtered_R1.fastq.trimmed.trimmed \
#    -o ${SM}_info_R1.fastq \
#    -p 1
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#rm ${SM}_filtered_R1.fastq.trimmed.trimmed
#
#${HOME}/utils/pairfq/Pairfq-0.17.0/bin/pairfq addinfo \
#    -i ${SM}_filtered_R2.fastq.trimmed.trimmed \
#    -o ${SM}_info_R2.fastq \
#    -p 2
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#rm ${SM}_filtered_R2.fastq.trimmed.trimmed
#
#${HOME}/utils/pairfq/Pairfq-0.14.3/bin/pairfq makepairs \
#    -f ${SM}_info_R1.fastq \
#    -r ${SM}_info_R2.fastq \
#    -fp ${SM}_sync_R1.fastq \
#    -rp ${SM}_sync_R2.fastq \
#    -fs ${SM}_unpaired_R1.fastq \
#    -rs ${SM}_unpaired_R2.fastq
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#rm ${SM}_info_R*.fastq ${SM}_unpaired_R*.fastq

#bwa mem -t ${BWANC} -M -R '@RG\tID:1\tLB:1\tPL:illumina\tSM:'${SM}'\tPU:1.1.1' \
#    ${BWAIX} \
#    <(zcat ${SM}_sync_R1.fastq.gz) \
#    <(zcat ${SM}_sync_R2.fastq.gz) | \
#    samtools sort -o ${SM}.raw.bam

#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

gatk MarkDuplicatesWithMateCigar \
    --INPUT ${SM}.raw.bam \
    --METRICS_FILE ${SM}.metrics.txt \
    --OUTPUT ${SM}.dedup.bam

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#rm ${SM}.raw.bam

gatk ReorderSam \
    --INPUT ${SM}.dedup.bam \
    --OUTPUT ${SM}.reorder.bam \
    --REFERENCE ${FASTA}

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.dedup.bam

samtools index ${SM}.reorder.bam

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

java -jar ${HOME}/utils/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
    --analysis_type RealignerTargetCreator \
    --reference_sequence ${FASTA} \
    --input_file ${SM}.reorder.bam \
    --known ${MILLS} \
    --known ${PHASE} \
    --out ${SM}.intervals \
    --filter_bases_not_stored \
    --intervals ${CPTUR}

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

java -jar ${HOME}/utils/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
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

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.reorder.bam ${SM}.reorder.bai ${SM}.intervals

gatk BaseRecalibrator \
    --input ${SM}.realign.bam \
    --known-sites ${DBSNP} \
    --known-sites ${MILLS} \
    --known-sites ${PHASE} \
    --output ${SM}.recal.table \
    --reference ${FASTA}

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

gatk ApplyBQSR \
    --bqsr-recal-file ${SM}.recal.table \
    --input ${SM}.realign.bam \
    --output ${SM}.recal.bam

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.realign.bam ${SM}.recal.table

samtools mpileup -AB \
    -f ${FASTA} \
    -o ${SM}.mpileup \
    ${SM}.recal.bam

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#rm ${SM}.recal.bam

java -jar ${HOME}/utils/varscan/VarScan-2.4.x/VarScan.v2.4.3.jar \
    mpileup2cns ${SM}.mpileup \
    --min-reads2 7 \
    --min-var-freq 0.0 \
    --p-value 0.05 \
    --output-vcf 1 \
    --variants \
    > ${SM}.vcf

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.mpileup

${HOME}/utils/annovar/convert2annovar.pl \
    --format vcf4 \
    --includeinfo \
    --outfile ${SM}.ann \
    ${SM}.vcf

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.vcf

${HOME}/utils/annovar/annotate_variation.pl \
    --outfile ${SM} \
    --buildver ${BUILD} \
    --hgvs \
    ${SM}.ann \
    ${ANNDB}

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.ann

${HOME}/src/Pipelines/variant_calling/format.py ${SM}.exonic_variant_function > ${SM}.xls

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.*variant_function ${SM}.log
exit 0
