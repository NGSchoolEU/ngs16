#!/bin/bash

#
# Copyright Russell S. Hamilton 2016
# rsh46@cam.ac.uk
#



#TMWG mkBS and oxBS performance
# >=99% for C conversion both oxBS and mkBS
# >=90% for hmC conversion in oxBS
# <10% for hmC conversion in mkBS (hmC/mkBS over­conversion error rate)
# <5% for mC conversion in both oxBS and mkBS (over­conversion error rate)


NUMREADS=1000000
ERRORRATE=0.25
BASEDIR="/Users/rhamilto/Documents/CTR-Teaching/NGSchool.eu/201608-MethylationCourse/"
DATADIR="NGSchool_GRCh38_Chr1_region"

# mkBS

Sherman \
  --length 100 \
  --number_of_seqs $NUMREADS \
  --genome_folder ${BASEDIR}/${DATADIR} \
  --paired_end \
  --minfrag 70 \
  --maxfrag 400 \
  --conversion_rate 99 \
  --error_rate ${ERRORRATE} \
  --variable_length_adapter 100


mv simulated_1.fastq mkbs_sim_${NUMREADS}_${ERRORRATE}_1.fastq; gzip mkbs_sim_${NUMREADS}_${ERRORRATE}_1.fastq
mv simulated_2.fastq mkbs_sim_${NUMREADS}_${ERRORRATE}_2.fastq; gzip mkbs_sim_${NUMREADS}_${ERRORRATE}_2.fastq

#oxBS

Sherman \
  --length 100 \
  --number_of_seqs $NUMREADS \
  --genome_folder ${BASEDIR}/${DATADIR} \
  --paired_end \
  --minfrag 70 \
  --maxfrag 400 \
  --conversion_rate 90 \
  --error_rate ${ERRORRATE} \
  --variable_length_adapter 100

mv simulated_1.fastq oxbs_sim_${NUMREADS}_${ERRORRATE}_1.fastq; gzip oxbs_sim_${NUMREADS}_${ERRORRATE}_1.fastq
mv simulated_2.fastq oxbs_sim_${NUMREADS}_${ERRORRATE}_2.fastq; gzip oxbs_sim_${NUMREADS}_${ERRORRATE}_2.fastq

