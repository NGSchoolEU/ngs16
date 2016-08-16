#!/bin/bash

peaks=`awk '{ SUM += $11 } END { print SUM }' $2`
total=`samtools view -c $1`

echo "$peaks / $total" | bc -l
