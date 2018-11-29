#!/usr/bin/env bash

ID=$1 # Input directory.
SF=$2 # Samples file (contains sample names).

OD=Results_UMI

mkdir -p ${OD}/FastQC

for SM in $(cat ${SF})
do
    I1=${ID}/${SM}_*_R1_001.fastq.gz
    I2=${ID}/${SM}_*_R2_001.fastq.gz
    SM=${OD}/${SM}

    /home/pcusco/utils/FastQC/fastqc ${I1} ${I2} --outdir Results_UMI/FastQC

    ./DS_part1_consensus.sh ${I1} ${I2} ${SM}

    ./DS_part2_map.sh ${SM}_read1_sscs.fq.gz ${SM}_read2_sscs.fq.gz ${SM}.sscs
    ./DS_part2_map.sh ${SM}_read1_dcs.fq.gz ${SM}_read2_dcs.fq.gz ${SM}.dcs

    ./DS_part3_sscs_vc.sh ${SM}.sscs 3
    ./DS_part3_dcs_vc.sh ${SM}.dcs

    ./DS_part4_annotate.sh ${SM}.sscs
    ./DS_part4_annotate.sh ${SM}.dcs
done

exit 0
