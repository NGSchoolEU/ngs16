#!/bin/bash

#
# Copyright Russell S. Hamilton 2016
# rsh46@cam.ac.uk
#

MODE=$1

picard_jar_location="/usr/local/picard-tools-2.2.4/picard.jar"
file=${MODE}"_sim_1000000_0.25_1_val_1_bismark_bt2_pe.srtd.bam"

java -jar $picard_jar_location CollectInsertSizeMetrics INPUT=$file OUTPUT=${file/.bam/}_picard_insert_size_metrics.txt HISTOGRAM_FILE=${file/.bam/}_picard_insert_size_plot.pdf METRIC_ACCUMULATION_LEVEL=ALL_READS

