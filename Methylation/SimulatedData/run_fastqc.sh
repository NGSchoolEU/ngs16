#!/bin/bash

#
# Copyright Russell S. Hamilton 2016
# rsh46@cam.ac.uk
#


MODE=$1

/Applications/FastQC.app/Contents/MacOS/fastqc -q ${MODE}_sim_1000000_0.25_1.fastq.gz 
/Applications/FastQC.app/Contents/MacOS/fastqc -q ${MODE}_sim_1000000_0.25_2.fastq.gz

