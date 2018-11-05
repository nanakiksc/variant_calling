#!/usr/bin/env bash

I1=$1
I2=$2
SM=$3

HOME=/home/pcusco
BWANC=8
BUILD=hg19
DATAD=${HOME}/data/${BUILD}
FASTA=${DATAD}/masked/${BUILD}.chr_sorted.fa
BWAIX=${DATAD}/masked/bwa/${BUILD}.masked
DBSNP=${DATAD}/dbsnp/All_20180423_chr.vcf.gz
MILLS=${DATAD}/gatk_resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
PHASE=${DATAD}/gatk_resources/1000G_phase1.indels.hg19.sites.vcf.gz
CPTUR=${HOME}/avivancos/plasma_umis/capture_regions/171106_HG19_VHIO_UMIs_EZ_HX3_capture_targets.bed
ANNDB=${HOME}/utils/annovar/humandb/

bwa aln -t ${BWANC} \
    ${BWAIX} \
    <(zcat ${I1}) \
    > ${SM}_R1.sai

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

bwa aln -t ${BWANC} \
    ${BWAIX} \
    <(zcat ${I2}) \
    > ${SM}_R2.sai

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

bwa sampe \
    -r '@RG\tID:RG1\tLB:DuplexSeq\tPL:illumina\tSM:'${SM}'\tPU:1.1.1' \
    ${BWAIX} \
    ${SM}_R1.sai \
    ${SM}_R2.sai \
    <(zcat ${I1}) \
    <(zcat ${I2}) | \
    samtools view -F 4 -q 1 -b | \
    samtools sort -o ${SM}.raw.bam

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}_R*.sai

gatk ReorderSam \
    --INPUT ${SM}.raw.bam \
    --OUTPUT ${SM}.reorder.bam \
    --REFERENCE ${FASTA}

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.raw.bam

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
rm ${SM}.reorder.bam* ${SM}.intervals

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
rm ${SM}.realign.ba* ${SM}.recal.table

samtools mpileup -B \
    -f ${FASTA} \
    -o ${SM}.mpileup \
    ${SM}.recal.bam

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
#rm ${SM}.recal.bam

java -jar ${HOME}/utils/varscan/VarScan-2.4.x/VarScan.v2.4.3.jar \
    mpileup2cns ${SM}.mpileup \
    --min-reads2 1 \
    --min-var-freq 0.0 \
    --p-value 0.05 \
    --strand-filter 0 \
    --output-vcf 1 \
    --variants \
    > ${SM}.vcf

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.mpileup

bgzip ${SM}.vcf

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

tabix ${SM}.vcf.gz

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

bcftools annotate \
    --annotations ${DBSNP} \
    --columns ID \
    --output ${SM}.dbsnp.vcf \
    --regions-file ${CPTUR} \
    ${SM}.vcf.gz

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.vcf.gz*

${HOME}/utils/annovar/convert2annovar.pl \
    --format vcf4 \
    --includeinfo \
    --outfile ${SM}.ann \
    ${SM}.dbsnp.vcf

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.dbsnp.vcf*

${HOME}/utils/annovar/annotate_variation.pl \
    --outfile ${SM} \
    --buildver ${BUILD} \
    --hgvs \
    ${SM}.ann \
    ${ANNDB}

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.ann

${HOME}/src/Pipelines/variant_calling/format.py ${SM} > ${SM}.xls

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.*variant_function ${SM}.log
exit 0
