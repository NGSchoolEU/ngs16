#!/bin/bash

#
# Copyright Russell S. Hamilton 2016
# rsh46@cam.ac.uk
#

MODE=$1

bismark_methylation_extractor --ignore_r2 1 --ignore_3prime_r2 2 --bedGraph --counts --gzip -p --no_overlap --report ${MODE}_sim_1000000_0.25_1_val_1_bismark_bt2_pe.deduplicated.bam

