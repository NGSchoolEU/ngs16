#!/bin/bash

#
# Copyright Russell S. Hamilton 2016
# rsh46@cam.ac.uk
#

MODE=$1


JAVA_OPTS="-Djava.awt.headless=true";
echo $JAVA_OPTS
qualimap bamqc -sd -c -bam ${MODE}_sim_1000000_0.25_1_val_1_bismark_bt2_pe.srtd.bam

