# NGSchool2016 materials
Materials prepared by the instructors of the [#NGSchool2016](https://ngschool.eu/). 

**Table of content**
- [Introduction](#Introduction)
- [Server setup](#server_setup)

## Introduction

In order to initialise all dependencies, execute:
```bash
source /ngschool/.bashrc
```

## Server setup
[Requested software](https://docs.google.com/spreadsheets/d/1uOQ2-1Yn_DyPd1_KFvpY-YmrS_87V6UZS19AG5akAQc/edit#gid=0)  
If you have no github account or no read permission for given repository, use `git clone https://` instead of `git clone git@`. 
```bash
# admin
sudo apt install git htop screen python-pip
sudo apt install libboost-iostreams-dev libboost-system-dev libboost-filesystem-dev zlib1g-dev libgsl0ldbl

sudo -H pip install -U pip

# biotools
sudo apt install fastqc soapdenovo2 ray velvet mummer bwa samtools bedtools igv fastx-toolkit last-align
sudo -H pip install -U numpy matplotlib biopython pysam

# tools for speakers
sudo apt install tophat bowtie bowtie2 idba ray trinity tabix picard-tools igv 

# first update R to 3.3+ https://www.r-bloggers.com/how-to-install-r-on-linux-ubuntu-16-04-xenial-xerus/
sudo apt install blast2 tigr-glimmer exonerate muscle fasttree mcl r-base

mkdir src && cd src
# bin
git clone git@github.com:lpryszcz/bin.git

# spades
wget http://spades.bioinf.spbau.ru/release3.9.0/SPAdes-3.9.0-Linux.tar.gz && tar xpfz SPAdes-3.9.0-Linux.tar.gz


# redundans
git clone git@github.com:lpryszcz/redundans.git
## install SSPACE & GapCloser
wget -q http://www.baseclear.com/base/download/41SSPACE-STANDARD-3.0_linux-x86_64.tar.gz
tar xpfz 41SSPACE-STANDARD-3.0_linux-x86_64.tar.gz
ln -s SSPACE-STANDARD-3.0_linux-x86_64 SSPACE
wget -q -O- http://cpansearch.perl.org/src/GBARR/perl5.005_03/lib/getopts.pl > SSPACE/dotlib/getopts.pl
wget -q http://downloads.sourceforge.net/project/soapdenovo2/GapCloser/bin/r6/GapCloser-bin-v1.12-r6.tgz
tar xpfz GapCloser-bin-v1.12-r6.tgz
rm 41SSPACE-STANDARD-3.0_linux-x86_64.tar.gz GapCloser-bin-v1.12-r6.tgz GapCloser_Manual.pdf
## test
source /ngschool/.bashrc
(cd redundans && ./redundans.py -v -i test/*.fq.gz -f test/contigs.fa -o test/run1 --sspacebin $SSPACEBIN)

# trinity & transdecoder
git clone git@github.com:trinityrnaseq/trinityrnaseq.git
(cd trinityrnaseq && make)
git clone git@github.com:TransDecoder/TransDecoder.git
(cd TransDecoder && make)

# cufflinks
wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar xpfz cufflinks-2.2.1.Linux_x86_64.tar.gz

# augustus
wget http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-3.0.2.tar.gz
tar xfpz augustus-3.0.2.tar.gz 
(cd augustus-3.0.2/src && make)

# GATK & snpEff

wget http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip && unzip snpEff_latest_core.zip

# Bismarc
git clone git@github.com:FelixKrueger/Bismark.git
wget http://smithlabresearch.org/downloads/preseq_linux_v2.0.tar.bz2 && tar -jxvf preseq_linux_v2.0.tar.bz2
wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.zip && unzip qualimap_v2.2.zip
wget http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip && unzip trim_galore_v0.4.1.zip

# 
sudo -H pip install git+https://github.com/ewels/MultiQC.git

git clone git@github.com:parseq/convector.git

# GATK https://software.broadinstitute.org/gatk/download/auth?package=GATK
(mkdir GATK && cd GATK && tar jxvf ../GenomeAnalysisTK-3.6.tar.bz2)

# InterProScan https://www.ebi.ac.uk/interpro/interproscan.html

# clean-up
rm *.zip *.gz *.bz2
```

```R
# jmarzec
# sudo R
source("http://bioconductor.org/biocLite.R"); biocLite("affy")
source("http://bioconductor.org/biocLite.R"); biocLite("affyPLM")
source("http://bioconductor.org/biocLite.R"); biocLite("arrayQualityMetrics")
source("http://bioconductor.org/biocLite.R"); biocLite("biomaRt")
source("http://bioconductor.org/biocLite.R"); biocLite("gcrma")
source("http://bioconductor.org/biocLite.R"); biocLite("geneplotter")
install.packages("gplots")
source("http://bioconductor.org/biocLite.R"); biocLite("hgu133acdf")
source("http://bioconductor.org/biocLite.R"); biocLite("hgu133plus2.db")
source("http://bioconductor.org/biocLite.R"); biocLite("limma")
source("http://bioconductor.org/biocLite.R"); biocLite("simpleaffy")
source("http://bioconductor.org/biocLite.R"); biocLite("sva")

# lpryszcz
source("http://bioconductor.org/biocLite.R"); biocLite("DESeq")
source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2")
source("http://bioconductor.org/biocLite.R"); biocLite("cummeRbund")

# rhamilton
install.packages("devtools")
library(devtools)
install_github("al2na/methylKit", build_vignettes=FALSE, 
  repos=BiocInstaller::biocinstallRepos(),
  dependencies=TRUE)

# mlapinski
install.packages('caTools')
source('http://bioconductor.org/biocLite.R'); biocLite('Rsamtools')
library(devtools); install_github("hms-dbmi/spp")

```
