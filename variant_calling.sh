#!/usr/bin/env bash

I1=$1
I2=$2
SM=$3

HOME=/home/pcusco
BWANC=8
VQUAL=30
BUILD=hg19
DATAD=${HOME}/data/${BUILD}
FASTA=${DATAD}/masked/${BUILD}.chr_sorted.fa
BWAIX=${DATAD}/masked/bwa/${BUILD}.masked
DBSNP=${DATAD}/dbsnp/All_20180423_chr.vcf.gz # chr
MILLS=${DATAD}/gatk_resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
PHASE=${DATAD}/gatk_resources/1000G_phase1.indels.hg19.sites.vcf.gz
CPTUR=${HOME}/avivancos/plasma_umis/capture_regions/171106_HG19_VHIO_UMIs_EZ_HX3_capture_targets.bed
ANNDB=${HOME}/utils/annovar/humandb/

fastq_quality_filter -q 30 -p 75 \
    -i <(zcat ${I1}) \
    -o ${SM}_R1.filtered.fastq

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

fastq_quality_filter -q 30 -p 75 \
    -i <(zcat ${I2}) \
    -o ${SM}_R2.filtered.fastq

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

${HOME}/utils/Homer/bin/homerTools trim \
    -5 TGACT \
    -minMatchLength 5 \
    -min 55 \
    ${SM}_R1.filtered.fastq \
    ${SM}_R2.filtered.fastq

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}_R*.filtered.fastq ${SM}_R*.filtered.fastq.lengths

${HOME}/utils/pairfq/Pairfq-0.17.0/bin/pairfq addinfo \
    -i ${SM}_R1.filtered.fastq.trimmed \
    -o ${SM}_R1.info.fastq \
    -p 1

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}_R1.filtered.fastq.trimmed

${HOME}/utils/pairfq/Pairfq-0.17.0/bin/pairfq addinfo \
    -i ${SM}_R2.filtered.fastq.trimmed \
    -o ${SM}_R2.info.fastq \
    -p 2

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}_R2.filtered.fastq.trimmed

${HOME}/utils/pairfq/Pairfq-0.14.3/bin/pairfq makepairs \
    -f ${SM}_R1.info.fastq \
    -r ${SM}_R2.info.fastq \
    -fp ${SM}_R1.sync.fastq \
    -rp ${SM}_R2.sync.fastq \
    -fs ${SM}_R1.unpaired.fastq \
    -rs ${SM}_R2.unpaired.fastq

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}_R*.info.fastq ${SM}_R*.unpaired.fastq

bwa aln -t ${BWANC} \
    ${BWAIX} \
    ${SM}_R1.sync.fastq \
    > ${SM}_R1.sai

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

bwa aln -t ${BWANC} \
    ${BWAIX} \
    ${SM}_R2.sync.fastq \
    > ${SM}_R2.sai

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

bwa sampe \
    -r '@RG\tID:RG1\tLB:DuplexSeq\tPL:illumina\tSM:'${SM}'\tPU:1.1.1' \
    ${BWAIX} \
    ${SM}_R1.sai \
    ${SM}_R2.sai \
    ${SM}_R1.sync.fastq \
    ${SM}_R2.sync.fastq | \
    samtools view -F 4 -q 1 -b | \
    samtools sort -o ${SM}.raw.bam

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}_R*.sync.fastq ${SM}_R*.sai

samtools index ${SM}.raw.bam

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

gatk MarkDuplicatesWithMateCigar \
    --INPUT ${SM}.raw.bam \
    --METRICS_FILE ${SM}.metrics.txt \
    --OUTPUT ${SM}.dedup.bam
#    --REMOVE_DUPLICATES true

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#rm ${SM}.raw.bam ${SM}.metrics.txt
rm ${SM}.metrics.txt

samtools index ${SM}.dedup.bam

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

#gatk ReorderSam \
#    --INPUT ${SM}.dedup.bam \
#    --OUTPUT ${SM}.reorder.bam \
#    --REFERENCE ${FASTA}
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
##rm ${SM}.dedup.bam
#
#samtools index ${SM}.reorder.bam
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#
#java -jar ${HOME}/utils/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
#    --analysis_type RealignerTargetCreator \
#    --reference_sequence ${FASTA} \
#    --input_file ${SM}.reorder.bam \
#    --known ${MILLS} \
#    --known ${PHASE} \
#    --out ${SM}.intervals \
#    --filter_bases_not_stored \
#    --intervals ${CPTUR}
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#
#java -jar ${HOME}/utils/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
#    --analysis_type IndelRealigner \
#    --reference_sequence ${FASTA} \
#    --input_file ${SM}.reorder.bam \
#    --targetIntervals ${SM}.intervals \
#    --out ${SM}.realign.bam \
#    --knownAlleles ${MILLS} \
#    --knownAlleles ${PHASE} \
#    --maxReadsForRealignment 150000 \
#    --intervals ${CPTUR} \
#    --filter_bases_not_stored
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
##rm ${SM}.reorder.bam* ${SM}.intervals
#rm ${SM}.intervals
#
#gatk BaseRecalibrator \
#    --input ${SM}.realign.bam \
#    --known-sites ${DBSNP} \
#    --known-sites ${MILLS} \
#    --known-sites ${PHASE} \
#    --output ${SM}.recal.table \
#    --reference ${FASTA}
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#
#gatk ApplyBQSR \
#    --bqsr-recal-file ${SM}.recal.table \
#    --input ${SM}.realign.bam \
#    --output ${SM}.recal.bam
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
##rm ${SM}.realign.ba* ${SM}.recal.table
#rm ${SM}.recal.table
#
#samtools mpileup -B \
#    -f ${FASTA} \
#    -o ${SM}.mpileup \
#    ${SM}.recal.bam
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
##rm ${SM}.recal.bam
#
#java -jar ${HOME}/utils/varscan/VarScan-2.4.x/VarScan.v2.4.3.jar \
#    mpileup2cns ${SM}.mpileup \
#    --min-reads2 7 \
#    --min-avg-qual 15 \
#    --min-var-freq 0.0 \
#    --p-value 0.05 \
#    --strand-filter 1 \
#    --output-vcf 1 \
#    --variants \
#    > ${SM}.vcf
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
##rm ${SM}.mpileup
#
#bgzip -f ${SM}.vcf
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#
#tabix ${SM}.vcf.gz
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#
#gatk VariantAnnotator \
#    --output ${SM}.dbsnp.vcf \
#    --variant ${SM}.vcf.gz \
#    --dbsnp ${DBSNP} \
#    --intervals ${CPTUR} \
#    --reference ${FASTA}
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
##rm ${SM}.vcf.gz*
#
#${HOME}/utils/annovar/convert2annovar.pl \
#    --format vcf4 \
#    --includeinfo \
#    --outfile ${SM}.ann \
#    ${SM}.dbsnp.vcf
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
##rm ${SM}.dbsnp.vcf*
#
#${HOME}/utils/annovar/annotate_variation.pl \
#    --outfile ${SM} \
#    --buildver ${BUILD} \
#    --hgvs \
#    ${SM}.ann \
#    ${ANNDB}
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
##rm ${SM}.ann
#
#${HOME}/src/Pipelines/variant_calling/format.py ${SM} > ${SM}.xls
#
#ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#rm ${SM}.*variant_function ${SM}.log

exit 0
