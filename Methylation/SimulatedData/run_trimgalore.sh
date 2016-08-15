#!/bin/bash

#
# Copyright Russell S. Hamilton 2016
# rsh46@cam.ac.uk
#


MODE=$1

trim_galore --paired --gzip ${MODE}_sim_1000000_0.25_1.fastq.gz ${MODE}_sim_1000000_0.25_2.fastq.gz

