# NGSchool2016 materials
Materials prepared by the instructors of the [#NGSchool2016](https://ngschool.eu/). 

Table of content
- [Server setup](#server_setup)


## Server setup
```bash
# admin
sudo apt install git htop screen python-pip
sudo -H pip install -U pip

# biotools
sudo apt install fastqc soapdenovo2 ray velvet mummer bwa samtools bedtools igv fastx-toolkit last-align
sudo -H pip install -U numpy matplotlib biopython pysam

# tools for speakers https://docs.google.com/spreadsheets/d/1uOQ2-1Yn_DyPd1_KFvpY-YmrS_87V6UZS19AG5akAQc/edit#gid=0
sudo apt install tophat bowtie bowtie2 idba ray trinity picard-tools igv blast2 tigr-glimmer exonerate muscle fasttree mcl r-base

# R for jmarzec / lpryszcz
sudo R
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
source("http://bioconductor.org/biocLite.R"); biocLite("DESeq")
source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2")


```
