#!/usr/bin/env bash

ID=$1 # Input directory.
SF=$2 # Samples file (contains sample names).

for SM in $(cat ${SF})
do
    I1=${ID}/${SM}_*_R1_001.fastq.gz
    I2=${ID}/${SM}_*_R2_001.fastq.gz

    mkdir -p FastQC
    /home/pcusco/utils/FastQC/fastqc ${I1} ${I2} --outdir FastQC

    ./DS_part1.sh ${I1} ${I2} ${SM}
    ./DS_part2.sh ${SM}_read1_sscs.fq.gz ${SM}_read2_sscs.fq.gz ${SM}.sscs
    ./DS_part2.sh ${SM}_read1_dcs.fq.gz ${SM}_read2_dcs.fq.gz ${SM}.dcs
done

exit 0
