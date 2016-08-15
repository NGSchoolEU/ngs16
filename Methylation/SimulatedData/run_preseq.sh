#!/bin/bash

#
# Copyright Russell S. Hamilton 2016
# rsh46@cam.ac.uk
#

MODE=$1

#preseq lc_extrap -B -P -e 4900100000 -o ${MODE}_sim_1000000_0.25_1_val_1_bismark_bt2_pe.srtd.preseq.lc_extrap.tsv ${MODE}_sim_1000000_0.25_1_val_1_bismark_bt2_pe.srtd.bam

preseq bound_pop -B -P -o ${MODE}_sim_1000000_0.25_1_val_1_bismark_bt2_pe.srtd.preseq.bound_pop.tsv ${MODE}_sim_1000000_0.25_1_val_1_bismark_bt2_pe.srtd.bam

