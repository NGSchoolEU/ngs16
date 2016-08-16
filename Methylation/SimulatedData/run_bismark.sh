#!/bin/bash

#
# Copyright Russell S. Hamilton 2016
# rsh46@cam.ac.uk
#

MODE=$1

bismark --bowtie2 ../NGSchool_GRCh38_Chr1_region/ -1 ${MODE}_sim_1000000_0.25_1_val_1.fq.gz -2 ${MODE}_sim_1000000_0.25_2_val_2.fq.gz

