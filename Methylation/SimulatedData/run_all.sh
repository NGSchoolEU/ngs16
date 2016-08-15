#!/bin/bash

#
# Copyright Russell S. Hamilton 2016
# rsh46@cam.ac.uk
#


#
# Genome preparation and make simulated reads
#

#./run_bismark_genome_prepare.sh
#./run_sherman.sh


#
# Process the bs-seq reads
#

MODE="mkbs"
./run_fastqc.sh ${MODE}
./run_trimgalore.sh ${MODE}
./run_bismark.sh ${MODE}
./run_bismark_dedup.sh ${MODE}
./run_bismark_extract.sh ${MODE}
./run_samtools_sort.sh ${MODE}
./run_samtools_index.sh ${MODE}
./run_preseq.sh ${MODE}
./run_qualimap.sh ${MODE}
./run_picard_is.sh ${MODE}

#
# Process the oxbs reads
#

MODE="oxbs"
./run_fastqc.sh ${MODE}
./run_trimgalore.sh ${MODE}
./run_bismark.sh ${MODE}
./run_bismark_dedup.sh ${MODE}
./run_bismark_extract.sh ${MODE}
./run_samtools_sort.sh ${MODE}
./run_samtools_index.sh ${MODE}
./run_preseq.sh ${MODE}
./run_qualimap.sh ${MODE}
./run_picard_is.sh ${MODE}

#
# Summarise the analysis
#

./run_multiqc.sh

#
# Compare bs and oxbs
#

Rscript NGSchool.eu.methylkit.R

