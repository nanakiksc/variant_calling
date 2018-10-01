#!/usr/bin/env bash

I1=$1 # Input file 1 (.fastq.gz).
I2=$2 # Input file 2 (.fastq.gz).
SM=$3 # Sample name.

mkdir -p FastQC
/home/pcusco/utils/FastQC/fastqc ${I1} ${I2} --outdir FastQC

./DS_part1.sh ${I1} ${I2} ${SM}
./DS_part2.sh ${SM}_read1_sscs.fq.gz ${SM}_read2_sscs.fq.gz ${SM}.sscs
./DS_part2.sh ${SM}_read1_dcs.fq.gz ${SM}_read2_dcs.fq.gz ${SM}.dcs
