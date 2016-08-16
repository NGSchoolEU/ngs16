#!/usr/bin/Rscript

#
# NGSchool.eu
#
# Copyright 2016 Russell S. Hamilton (rsh46@cam.ac.uk)


#
# To install Methyl=kit uncomment the section below
#
#library(devtools)
# Standard Release
#install_github("al2na/methylKit", build_vignettes=FALSE, repos=BiocInstaller::biocinstallRepos(), dependencies=TRUE)
# Development Version
#install_github("al2na/methylKit", build_vignettes=FALSE, repos=BiocInstaller::biocinstallRepos(), ref="development", dependencies=TRUE)




library(methylKit)

#
# Read in the bismark generated coverage files
# readBismarkCoverage and fread code taken from: https://gist.github.com/al2na/4839e615e2401d73fe51


readBismarkCoverage<-function( location,sample.id,assembly="unknown",treatment, context="CpG",min.cov=10)
#' Read bismark coverage file as a methylKit object
#' 
#' Bismark aligner can output methylation information per base in
#' multiple different formats. This function reads coverage files,
#' which have chr,start,end, number of cytosines (methylated bases) 
#' and number of thymines (unmethylated bases).
#' 
#' @param location a list or vector of file paths to coverage files
#'     
#' @param sample.id a list or vector of sample ids
#' @param assembly a string for genome assembly. Any string would work.
#' @param treatment if there are multiple files to be read a treatment 
#'                  vector should be supplied.
#' @param context a string for context of methylation such as: "CpG" or "CHG"
#' @param min.cov a numeric value for minimum coverage. Bases that have coverage
#' below this value will be removed.
#' 
#' @return methylRaw or methylRawList objects
{
  if(length(location)>1)
    {
      stopifnot(length(location)==length(sample.id),length(location)==length(treatment))
    }
  
  result=list()
  for(i in 1:length(location))
     {
       df=fread.gzipped(location[[i]],data.table=FALSE)
       # remove low coverage stuff
       df=df[ (df[,5]+df[,6]) >= min.cov ,]
       # make the object (arrange columns of df), put it in a list
       result[[i]]= new("methylRaw",data.frame(chr=df[,1],start=df[,2],end=df[,3],strand="*",coverage=(df[,5]+df[,6]),numCs=df[,5],numTs=df[,6]),
                       sample.id=sample.id[[i]], assembly=assembly,context=context,resolution="base")
     }
  
  if(length(result) == 1)
    {
      return(result[[1]])
    }
  else
    {
      new("methylRawList",result,treatment=treatment)
    }
}

# reads gzipped files,
fread.gzipped<-function(filepath,...){
  require(R.utils)
  require(data.table)
  if (R.utils::isGzipped(filepath))
    {
      if(.Platform$OS.type == "unix") 
        {
          filepath=paste("gunzip -c",filepath)
        } 
      else 
        {
          filepath <- R.utils::gunzip(filepath,temporary = FALSE, overwrite = TRUE, remove = FALSE)
        }
    }
  fread(filepath,...)
}


print("Read in the bismark generated coverage files:")
bismark_cov_files=list("mkbs_sim_1000000_0.25_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz", "oxbs_sim_1000000_0.25_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")

print("Creating a methyl-kit object")
methylkit.obj = readBismarkCoverage(bismark_cov_files,  sample.id=list( "mkbs","oxbs"), assembly="hg38", treatment=c(1,0), min.cov=5)

print("Creating a united object and sample clustering")
meth <- unite(methylkit.obj)
hc   <- clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)

print("Creating a ward correlation plot")
pdf(file="NGSchool.eu.methylkit.PCASamples.ward_corr_plot.pdf");
clusterSamples(meth,dist="correlation", method="ward", plot=TRUE)
dev.off()

print("Creating a Correlation Plot")
pdf(file="NGSchool.eu.methylkit.CorrelationPlot.pdf");
getCorrelation(meth,plot=T)
dev.off()

print("Creating a PCA Scree Plot")
pdf(file="NGSchool.eu.methylkit.PCASamples.screeplot.pdf")
PCASamples(meth,screeplot=TRUE)
dev.off()

print("Creating a PCA Plot")
pdf(file="NGSchool.eu.methylkit.PCASamples.pdf")
PCASamples(meth)
dev.off()


print("Calculating differential metylation")

diffMeth=calculateDiffMeth(meth)
write.table(diffMeth, file="NGSchool.eu.methylkit.DiffMeth.tsv", sep='\t', quote=FALSE)

diffMeth.hyper = getMethylDiff(diffMeth,differenc=1,qvalue=0.1,type="hyper")
write.table(diffMeth.hyper,"NGSchool.eu.methylkit.hyper_methylated.tsv",sep='\t', quote=FALSE)

diffMeth.hypo = getMethylDiff(diffMeth,differenc=1,qvalue=0.1,type="hypo")
write.table(diffMeth.hypo,"NGSchool.eu.methylkit.hypo_methylated.tsv",sep='\t', quote=FALSE)

diffMeth = getMethylDiff(diffMeth,differenc=1,qvalue=0.1)
write.table(diffMeth,"NGSchool.eu.methylkit.differentialy_methylated.tsv",sep='\t', quote=FALSE)

pdf("NGSchool.eu.methylkit.diffMethPerChr.pdf");
diffMethPerChr(diffMeth,plot=TRUE,qvalue.cutoff=0.1,meth.cutoff=1)
dev.off()


#--------------------------------------------------------------------------------------------------
# END
#--------------------------------------------------------------------------------------------------
