#!/usr/bin/env bash

I1=$1
I2=$2
SM=$3

HOME=/home/pcusco
TAGLN=12 # 6
REPFL=9 # 5
BUILD=hg19

fastq_quality_filter -q 30 -p 75 \
    -i <(zcat ${I1}) \
    -o ${SM}_R1.filtered.fastq

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

fastq_quality_filter -q 30 -p 75 \
    -i <(zcat ${I2}) \
    -o ${SM}_R2.filtered.fastq

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi

${HOME}/utils/Homer/bin/homerTools trim \
    -3 AGATCGGAAGAGCACACGTCT \
    -mis 2 \
    -minMatchLength 4 \
    -min 55 \
    ${SM}_R1.filtered.fastq \
    ${SM}_R2.filtered.fastq

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}_R*.filtered.fastq ${SM}_R*.filtered.fastq.lengths

${HOME}/utils/pairfq/Pairfq-0.17.0/bin/pairfq addinfo \
    -i ${SM}_R1.filtered.fastq.trimmed\
    -o ${SM}_R1.info.fastq \
    -p 1

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}_R1.filtered.fastq.trimmed

${HOME}/utils/pairfq/Pairfq-0.17.0/bin/pairfq addinfo \
    -i ${SM}_R2.filtered.fastq.trimmed\
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

gatk FastqToSam \
    --FASTQ ${SM}_R1.sync.fastq \
    --FASTQ2 ${SM}_R2.sync.fastq \
    --OUTPUT ${SM}.unmapped.bam \
    --SAMPLE_NAME ${SM}

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}_R*.sync.fastq

python ${HOME}/src/Python/Duplex-Sequencing/UnifiedConsensusMaker.py \
    --input ${SM}.unmapped.bam \
    --taglen ${TAGLN} \
    --tagstats \
    --write-sscs \
    --rep_filt ${REPFL} \
    --prefix ${SM}

ec=$?; if [ $ec -ne 0 ]; then exit $ec; fi
rm ${SM}.unmapped.bam
