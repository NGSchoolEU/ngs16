#!/bin/bash

#
# Copyright Russell S. Hamilton 2016
# rsh46@cam.ac.uk
#

MODE=$1

deduplicate_bismark -p --bam ${MODE}_sim_1000000_0.25_1_val_1_bismark_bt2_pe.bam
