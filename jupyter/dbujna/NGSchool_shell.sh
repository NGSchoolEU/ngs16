#!/bin/bash

# Change to the chip_seq directory and type:
source ./.chipseqrc

# FastQC analysis
# To asses the quality of the sequencing experiment run:
mkdir fastqc
fastqc -t 2 -o "./fastqc" reads/*.fastq.gz

# Cutadapt adapter trimming
#
# To trim the adapters and low quality bases we can use Cutadapt:
# http://cutadapt.readthedocs.io/en/stable/guide.html
#
# The following options are specified:
# -f ; input file format
# -e ; error rate
# -q ; quality cutoff
# -O ; minimum adapter overlap
# -a ; 3' adapter sequence
# -m ; minimum length
# -o ; output file [STDOUT]

# To extract the name of our file we will use Bash parameter substitution
# (For further reading and other tricks refer to 
# http://tldp.org/LDP/abs/html/parameter-substitution.html#PARAMSUBREF )

for f in reads/*.fastq.gz; do
  name=${f%.*.*}; name=${name##*/};
  out_dir=cutadapt
  out_file=$out_dir/trimmed/$name.trimmed.fq.gz
  report_file=$out_dir/$name.report.txt
  if [ ! -d $out_dir/trimmed ]; then mkdir -v -p $out_dir/trimmed; fi
  if [ ! -d fastqc/trimmed ]; then mkdir -v -p fastqc/trimmed; fi
  date +'%D %T'; echo $name;
  cutadapt \
    -f fastq \
    -e 0.1 \
    -q 20 \
    -O 3 \
    -a AGATCGGAAGAGC \
    -m 20 \
    -o "$out_file" \
    $f | tee $report_file;
  fastqc -o "fastqc/trimmed" $out_file
done

# Align the reads

# prepare the index
mkdir -pv ./ref/danRer10.bwa.index/
bwa index -p "./ref/danRer10.bwa.index/danRer10.bwa" ./ref/danRer10.dna_sm.toplevel.fa

# Do the proper alignment now

# To obtain high quality reads for further processing we:
# 1. First align the reads to the reference genome with the bwa mem tool
# 2. Next we use the shell pipe to redirect the resulting SAM alignment to the samtools tool with wich we filter out the low quality alignments 
# 3. Subsequently we redirect the output of the previous command and sort the good quality alignments to store it in a BAM file
# 4. In the end we index the sorted BAM file for fast access
# 
# With samtools we will filter out the reads having the SAM flag value:
# 4 - read unmapped
# 256 - read is not primary alignment
# For processing we can add those flag values to each other.

for f in ./reads/*.fastq.gz; do
  outDir=./bwa;
  name=`basename $f`; name=${name%%.*};
  if [[ ! -d "$outDir" ]]; then mkdir -vp $outDir; fi
  if [[ ! -s "$outDir/$name.sorted.bam" ]]; then
    echo; date +'%D %T'; echo $name; echo;
    bwa mem \
      -t 2 \
      -M \
      ./ref/danRer10.bwa.index/danRer10.bwa \
      $f \
      | samtools view -b -u -q10 -F260 - \
      | samtools sort -o $outDir/$name.sorted.bam -O bam -T $outDir/$name.tmp -@ 2 -;
    samtools index $outDir/$name.sorted.bam;
  fi
done
# for the old version of samtools in the Ubuntu repository, substitute the last 2 commands with:
#| samtools view -S -b -h -q10 -F260 - \
#| samtools sort -f - $outdir/$name.sorted.bam;

# How many high quality ("unique") alignments we will wor with:
for f in ./bwa/*.bam; do echo $f; samtools view -c $f; done

# Calculate the NRF 
# We can simply calculate the Non-redundant Read Fraction using samtools and UNIX shell:

# specify the file
file=./bwa/input.sorted.bam
# and calculate the read counts
total=`samtools view -c $file`;
# We need to extract the coordinates of the reads aligning to the Watson and Crick strand
# separately (do you know why?). To do this we will check for the occurence of the flag value 16.
plusStrand=`samtools view -F16 $file | cut -f3-4 | sort -u | wc -l`;
minusStrand=`samtools view -f16 $file | cut -f3-4 | sort -u | wc -l`;

NRF=`echo "($plusStrand + $minusStrand)/$total" | bc -l`

echo $NRF

# Be aware that the NRF value will be lower the higher the sequencing depth.
# We could run this for all the files or do it in python in a bit more sophisticated way, with random read sampling.

# NSC and RSC calculation
#
# To check the cross-correlation of the reads aligned to both strands we will run the SPP.
mkdir -v ./spp;
for f in ./bwa/*.bam; do
  Rscript ./.local/share/phantompeakqualtools/run_spp_modified.R \
    -c=$f \
    -s=-100:2:600 \
    -savp \
    -odir=./spp \
    -out=./spp/spp.result
done

# Signal extraction scaling

# To perform the first two steps of the signal extraction scaling analysis we will use bedtools programs.
# It has a great documentation under:
# http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html
# First we will split our reference genome into 1kb non-overlapping windows:
mkdir bedtools;
bedtools makewindows -g ./ref/danRer10.83.len -w 1000 >./bedtools/windows.bed;
# Next we can calculate the read coverage of each of this window. We will use the following command:
for f in ./bwa/*.sorted.bam; do
  name=`basename $f`; name=${name%%.*};
  bedtools coverage \
    -a ./bedtools/windows.bed \
    -b $f \
    -counts \
    -sorted \
    -g ./ref/danRer10.83.len \
    >./bedtools/${name}_coverage.bed
done
  
# Now we can proceed with the next steps in python!

# ---

# Call the peaks with MACS2 peak caller

mkdir ./macs2;
macs2 callpeak \
  --treatment "./bwa/sox2_chip.sorted.bam" \
  --control "./bwa/input.sorted.bam" \
  --name "sox2_chip" \
  --outdir "./macs2" \
  --format BAM \
  --gsize 1.01e8 \
  --qvalue 0.01 \
  --bdg \
  2>&1 \
  | tee ./macs2/sox2_chip.macs2.log

# Create the MACS2 diagnostic plots
cd ./macs2;
Rscript sox2_chip_model.r;
cd ..;

# To see the signal enrichment in the genome browser we need to create files in the BigWig format. To do this we first need to compare
# the signal in the ChIP sample to the signal in the input sample. We will do it in MACS2, using fold enrichment as the output method.
macs2 bdgcmp -t ./macs2/sox2_chip_treat_pileup.bdg -c ./macs2/sox2_chip_control_lambda.bdg -o ./macs2/sox2_chip_FE.bdg -m FE
# Next we will change the Bedgraph files to the BigWig files needed for viewing with the following command:
bdg2bw ./macs2/sox2_chip_FE.bdg ./ref/danRer10.83.len

# Now we can run the IGV, load the provided igv.genome file (in the .ref/ folder), our alignments and predicted peaks to see the results of our analysis.

# Calculate FRiP

# To to calculate the fraction of reads in peaks we can again use bedtools.

bedtools coverage \
  -a ./macs2/sox2_chip_peaks.narrowPeak \
  -b ./bwa/sox2_chip.sorted.bam \
  -counts \
  -sorted \
  -g ./ref/danRer10.83.len \
  >./bedtools/sox2_chip_peak_coverage.bed

# We can sum all the read counts from the 11th column using awk, for simplicity just run the provided script:
calculateFRiP.sh ./bwa/sox2_chip.sorted.bam ./bedtools/sox2_chip_peak_coverage.bed

# In the end we can obtain the list of genes closest to our detected peaks:
bedtools closest \
  -a ./macs2/sox2_chip_peaks.narrowPeak \
  -b ./ref/ensembl_genes.sorted.bed \
  >./macs2/sox2_chip.genes.closest;

# vim: tw=0
