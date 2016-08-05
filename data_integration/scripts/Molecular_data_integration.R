################################################################################
#
#   File name: Molecular_data_integration.R
#
#   Authors: Jacek Marzec ( j.marzec@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
################################################################################

################################################################################
#
#	Description: Script for integrating prostate cancer expression data generated with Affymetrix U133 Plus 2.0 array
#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

#===============================================================================
#    Functions
#===============================================================================

##### Assign colours to analysed datasets
getDatasetsColours <- function(datasets) {
    
    ##### Predefined selection of colours for datasets
    datasets.colours <- c("bisque","orange","firebrick","lightslategrey","darkseagreen","darkcyan","dodgerblue")
    
    f.datasets <- factor(datasets)
    vec.datasets <- datasets.colours[1:length(levels(f.datasets))]
    datasets.colour <- rep(0,length(f.datasets))
    for(i in 1:length(f.datasets))
    datasets.colour[i] <- vec.datasets[ f.datasets[i] == levels(f.datasets)]
    
    return( list(vec.datasets, datasets.colour) )
}


##### Assign colours to analysed groups
getTargetsColours <- function(targets) {
    
    ##### Predefined selection of colours for groups
    targets.colours <- c("red","blue","green","darkgoldenrod","darkred","deepskyblue")
    
    f.targets <- factor(targets)
    vec.targets <- targets.colours[1:length(levels(f.targets))]
    targets.colour <- rep(0,length(f.targets))
    for(i in 1:length(f.targets))
    targets.colour[i] <- vec.targets[ f.targets[i] == levels(f.targets)]
    
    return( list(vec.targets, targets.colour) )
}


################################################################################
#    Load libraries
################################################################################

library(affy)
library(sva)
library(arrayQualityMetrics)
library(gcrma)
library(limma)
library(hgu133plus2.db)
library(gplots)


##### Source script for coloured dendrogram
source("/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/scripts/a2R_code.R")


################################################################################
#
#    QUALITY CONTROL:
#
#    Load data and perform QC for each dataset separately
#
################################################################################

#===============================================================================
#    Load data from Affymetrix 'CEL' files for dataset 1 (GSE17951)
#===============================================================================

##### Change working directory to data folder
setwd("/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data")


##### Read annotation file
pd <- read.AnnotatedDataFrame(filename = "target_GSE17951.txt" ,header = TRUE, row.name = "Name", sep = "\t")


##### Read CEL files into an Affybatch
dat <- ReadAffy(filenames = pd$FileName, sampleNames = sampleNames(pd), phenoData = pd)


#===============================================================================
#    Data QC
#===============================================================================

##### Data quality control using arrayQualityMetrics package
##### NOTE: Computationally intense step, may take few minutes to complete
arrayQuality <- arrayQualityMetrics(expressionset = dat, outdir = "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/QC/GSE17951", reporttitle = "arrayQualityMetrics report for GSE17951", force = TRUE, do.logtransform = TRUE, intgroup = "Target")


##### Report detected samples to excluded from downstream analysis
samples2exclude <- NULL


#===============================================================================
#    Load data from Affymetrix 'CEL' files for dataset 2 (GSE55945)
#===============================================================================

##### Read annotation file
pd <- read.AnnotatedDataFrame(filename = "target_GSE55945.txt", header = TRUE, row.name = "Name", sep = "\t")


##### Read CEL files into an Affybatch
dat <- ReadAffy(filenames = pd$FileName, sampleNames = sampleNames(pd), phenoData = pd)


#===============================================================================
#    Data QC
#===============================================================================

##### Data quality control using arrayQualityMetrics package
##### NOTE: Computationally intense step, may take few minutes to complete
arrayQuality <- arrayQualityMetrics(expressionset = dat, outdir = "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/QC/GSE55945", reporttitle = "arrayQualityMetrics report for GSE55945", force = TRUE, do.logtransform = TRUE, intgroup = "Target")


##### Report detected samples to excluded from downstream analysis
samples2exclude <- c(samples2exclude, "GSE55945_3", "GSE55945_4", "GSE55945_6", "GSE55945_15")


#===============================================================================
#    Load data from Affymetrix 'CEL' files for dataset 3 (GSE3325)
#===============================================================================

##### Read annotation file
pd <- read.AnnotatedDataFrame(filename = "target_GSE3325.txt", header = TRUE, row.name = "Name", sep = "\t")


##### Read CEL files into an Affybatch
dat <- ReadAffy(filenames = pd$FileName, sampleNames = sampleNames(pd), phenoData = pd)


#===============================================================================
#    Data QC
#===============================================================================

##### Data quality control using arrayQualityMetrics package
##### NOTE: Computationally intense step, may take few minutes to complete
arrayQuality <- arrayQualityMetrics(expressionset = dat, outdir = "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/QC/GSE3325", reporttitle = "arrayQualityMetrics report for GSE3325", force = TRUE, do.logtransform = TRUE, intgroup = "Target")


##### Report detected samples to excluded from downstream analysis
samples2exclude <- c(samples2exclude, "GSE3325_9")




#===============================================================================
#    Automated outliers detection (additional step)
#===============================================================================

##### Compute useful summary statistics from a data object
preparedData <- prepdata(expressionset = dat, intgroup = "Target", do.logtransform = TRUE)


##### Array intensity distributions
QCboxplot <- aqm.boxplot(preparedData)
QCboxplot@outliers


##### Between array comparison
QCheatmap <- aqm.heatmap(preparedData)
QCheatmap@outliers


##### MA plots
QCmaplot <- aqm.maplot(preparedData)
QCmaplot@outliers



##### For Affymetrix specific sections
preparedAffy <- prepaffy(expressionset = dat, x = preparedData)


##### RLE
QCrle <- aqm.rle(preparedAffy,)
QCrle@outliers


##### NUSE
QCnuse <- aqm.nuse(preparedAffy)
QCnuse@outliers


##### Customising arrayQualityMetrics report
qm <- list("Boxplot" = QCboxplot, "Heatmap" = QCheatmap, "MAplot" = QCmaplot, "RLE" = QCrle, "NUSE" = QCnuse)

aqm.writereport(modules = qm, reporttitle = "arrayQualityMetrics report for GSE3325", outdir = "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/QC/GSE3325_custom", arrayTable = pData(dat))




################################################################################
#
#    PRE-PROCESSING (1):
#
#    Normalisation of all datasets collectively using GC-RMA
#
################################################################################

#===============================================================================
#     Assign colours to datasets and groups
#===============================================================================

##### Read file with information for all datasets
datasetsInfo <- read.table(file = "target_all.txt", sep = "\t", as.is = TRUE, header = TRUE, row.names = 1)

##### Ignore samples identified as outliers
datasetsInfo <- datasetsInfo[setdiff(rownames(datasetsInfo), as.vector(samples2exclude)),]


##### Assign different colours for samples from individual datasets
datasets <- datasetsInfo$Dataset

datasets.colour <- getDatasetsColours(datasets)


##### Assign different colours for samples representing individual groups
targets <- datasetsInfo$Target

targets.colour <- getTargetsColours(targets)


##### Record number of datasets
datasets_No <- max(as.numeric(factor(datasets)))


#===============================================================================
#    Load data from Affymetrix 'CEL' files for all datasets
#===============================================================================

##### Read annotation files for all datasets
pd <- read.AnnotatedDataFrame(filename = "target_all.txt", header = TRUE, row.name = "Name", sep = "\t")

##### Ignore samples identified as outliers
pd <- pd[setdiff(sampleNames(pd), as.vector(samples2exclude))]


##### Read CEL files into an Affybatch
dat <- ReadAffy(filenames = pd$FileName, sampleNames = sampleNames(pd), phenoData = pd)



##### Create folder for results
system("mkdir /Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/results")

##### Change working directory to results folder
setwd("/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/results")



#===============================================================================
#     Normalisation with GC-RMA
#===============================================================================

##### Perform normalisation with GC-RMA method
##### NOTE: Computationally intensive step, it may take few minutes to complete
gcrma = gcrma(dat)


##### Generate signal intensity before and after GC-RMA normalisation
pdf("Norm_density.pdf")

#####  Before normalisation
hist(dat, col = datasets.colour[[2]], lty = 1)

#####  After normalisation
hist(gcrma, col = datasets.colour[[2]], lty = 1)
dev.off()


##### Generate box-and-whisker plot
pdf("Norm_boxplot.pdf", pointsize = 8, width = 0.2*length(pd$FileName), height = 6)
par(mar=c(13, 4, 3, 2))

#####  Before normalisation
boxplot(dat, col = datasets.colour[[2]], las = 2)

#####  After normalisation
boxplot(exprs(gcrma), col = datasets.colour[[2]], main = "GC-RMA normalised", las = 2)
dev.off()


##### Write normalised expression data into a file
write.exprs(gcrma, file = "Norm_data.exp", sep = "\t")


################################################################################
#
#    PRE-PROCESSING (2):
#
#    Study effects assessment
#
################################################################################

#===============================================================================
#     Unsupervised clustering
#===============================================================================

#####  Extract matrix of expression values
data <- exprs(gcrma)


#####  Compute distance matrix and hierarchical clustering
d.usa <- dist(t(data), method = "euc")
h.usa <- hclust(d.usa, method = "ward.D2")


#####  Generate coloured dentrogram (indicate datasets)
pdf("Study_effects_cluster_datasets.pdf", width = 0.3*ncol(data), height = 6)
h.usa$labels = colnames(data)
par(mar = c(2,2,2,6))
A2Rplot(h.usa, fact.sup = datasets, box = FALSE, show.labels = TRUE, col.up = "black", main = "", k = length(levels(factor(datasets))) ) # k = changes the detail/number of subgroups shown
dev.off()


#####  Generate coloured dentrogram (indicate groups)
pdf("Study_effects_cluster_targets.pdf", width = 0.3*ncol(data), height = 6)
h.usa$labels = colnames(data)
par(mar = c(2,2,2,6))
A2Rplot(h.usa, fact.sup = targets, box = FALSE, show.labels=TRUE, col.up = "black", main="", k = length(levels(factor(targets))) ) # k = changes the detail/number of subgroups shown.
dev.off()


#####  Generate dentrogram
pdf("Study_effects_dendrogram.pdf", width = 0.2*ncol(data)+2, height = 6, pointsize = 12)
par(mar = c(2,5,2,0))
plot(h.usa, xlab = "", labels = paste(colnames(data), targets, sep = "       "), hang = -1, main="")
dev.off()



#===============================================================================
#     Principal components analysis
#===============================================================================

#####  Keep only probes with variance > 0 across all samples
rsd <- apply(data, 1, sd)
data <- data[rsd > 0,]


#####  Perform PCA
data_pca <- prcomp(t(data), scale=TRUE)


#####  Get variance importance for all principal components
importance_pca <- summary(data_pca)$importance[2,]
importance_pca <- paste(round(100*importance_pca, 2), "%", sep="")


#####  Set point symbols so that each represents single dataset on the PCA plots
pchs <- rep(list(c(16,1),c(15,0),c(17,2),c(16,10),c(15,7),c(16,13),c(16,1)),10)


#####  Generate PCA plot
pdf("Study_effects_PCA.pdf")
plot(data_pca$x[,1], data_pca$x[,2], type = "n", xlab = paste("PC1 (",importance_pca[1],")", sep = ""), ylab = paste("PC2 (",importance_pca[2],")", sep = ""), ylim =c ( min(data_pca$x[,2]),max(data_pca$x[,2]) + (abs(min(data_pca$x[,2])) + abs(max(data_pca$x[,2])))/4 ))

#####  Use different shape for each dataset
for (i in 1:datasets_No) {
    points(data_pca$x[,1][as.numeric(factor(datasets)) == i], data_pca$x[,2][as.numeric(factor(datasets)) == i], pch = pchs[[i]][1], col = targets.colour[[2]][as.numeric(factor(datasets)) == i])
    points(data_pca$x[,1][as.numeric(factor(datasets)) == i], data_pca$x[,2][as.numeric(factor(datasets)) == i], pch = pchs[[i]][2], col = "black")
}

#####  Adding the legend
legend("topleft", legend = levels(factor(targets)), pch = 16, col = targets.colour[[1]], box.col = "transparent")
legend("topleft", legend = c(rep("", length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
dev.off()



################################################################################
#
#    INTEGRATIVE ANALYSIS (1):
#
#    Data adjustment for study effects
#
################################################################################


#####  Specify known batches
batch <-as.vector(datasets)

#####  Create model matrix for outcome of interest and other covariates beside batch
f.model <- factor(targets, levels = levels(factor(targets)))

mod <- model.matrix(~f.model)


#####  Adjust data for batch effects using ComBat
data_combat <- ComBat(dat = data, batch = batch, mod = mod)


#===============================================================================
#     Boxplot of combined datasets
#===============================================================================

pdf("Study_effects_boxplot.pdf", pointsize = 8, width = 0.2*ncol(data), height = 6)
par(mar = c(13, 4, 3, 2))

#####  Before adjusting for study effects
boxplot(data.frame(data), col = datasets.colour[[2]], las = 2)

#####  After adjusting for study effects
boxplot(data.frame(data_combat), col = datasets.colour[[2]], main = "Batch effect adjusted", las = 2)
dev.off()


#===============================================================================
#     Unsupervised clustering
#===============================================================================

#####  Compute distance matrix and hierarchical clustering
d.usa <- dist(t(data_combat), method = "euc")
h.usa <- hclust(d.usa, method = "ward.D2")


##### Generate coloured dentrogram (indicate datasets)
pdf("Study_effects_cluster_datasets_ComBat.pdf", width = 0.3*ncol(data), height = 6)
h.usa$labels = colnames(data_combat)
par(mar = c(2,2,2,6))
A2Rplot(h.usa, fact.sup = datasets, box = FALSE, show.labels = TRUE, col.up = "black", main="", k=length(levels(factor(datasets))) ) # k = changes the detail/number of subgroups shown.
dev.off()


##### Generate coloured dentrogram (indicate groups)
pdf("Study_effects_cluster_targets_ComBat.pdf", width = 0.3*ncol(data), height = 6)
h.usa$labels = colnames(data_combat)
par(mar = c(2,2,2,6))
A2Rplot(h.usa, fact.sup = targets, box = FALSE, show.labels = TRUE, col.up = "black", main=" ", k = length(levels(factor(targets))) ) # k = changes the detail/number of subgroups shown.
dev.off()


##### Generate dentrogram
pdf(paste("Study_effects_dendrogram_ComBat.pdf", sep = ""), width = 0.2*ncol(data)+2, height = 6, pointsize = 12)
par(mar = c(2,5,2,0))
plot(h.usa, xlab = "", labels = paste(colnames(data_combat), targets, sep="       "), hang = -1, main = "")
dev.off()



#===============================================================================
#     Principal components analysis
#===============================================================================

#####  Keep only probes with variance > 0 across all samples
rsd <- apply(data_combat, 1, sd)
data_combat <- data_combat[rsd > 0,]


#####  Perform PCA
data_combat_pca <- prcomp(t(data_combat), scale=TRUE)


#####  Get variance importance for all principal components
importance_pca <- summary(data_combat_pca)$importance[2,]
importance_pca <- paste(round(100*importance_pca, 2), "%", sep="")


#####  Generate PCA plot
pdf("Study_effects_PCA_ComBat.pdf")
plot(data_combat_pca$x[,1], data_combat_pca$x[,2], type = "n", xlab = paste("PC1 (",importance_pca[1],")", sep = ""), ylab = paste("PC2 (",importance_pca[2],")", sep = ""), ylim = c( min(data_combat_pca$x[,2]), max(data_combat_pca$x[,2]) + (abs(min(data_combat_pca$x[,2]))+ abs(max(data_combat_pca$x[,2])))/4 ))

#####  Use different shape for each dataset
for (i in 1:datasets_No) {
    points(data_combat_pca$x[,1][as.numeric(factor(datasets)) == i], data_combat_pca$x[,2][as.numeric(factor(datasets)) ==i ], pch = pchs[[i]][1], col = targets.colour[[2]][as.numeric(factor(datasets)) == i])
    points(data_combat_pca$x[,1][as.numeric(factor(datasets)) == i], data_combat_pca$x[,2][as.numeric(factor(datasets)) == i], pch = pchs[[i]][2], col = "black")
}

#####  Adding the legend
legend("topleft", legend = levels(factor(targets)), pch = 16, col = targets.colour[[1]], box.col = "transparent")
legend("topleft", legend = c(rep("", length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
dev.off()



################################################################################
#
#    INTEGRATIVE ANALYSIS (2):
#
#    Differential expression analysis
#
################################################################################

#===============================================================================
#     Non-specific filtering
#===============================================================================

#####  Use 60% of all genes with the highest expression variance across samples in the non-specific filtering step
filterThreshold <- round(0.6*nrow(data_combat), digits = 0)


rsd <- apply(data_combat, 1, sd)
sel <- order(rsd, decreasing=TRUE)[1:filterThreshold]
data_combat <- data_combat[sel,]


#===============================================================================
#     Differenatial expression analysis
#===============================================================================

#####  Create a model matrix representing the study design, with rows corresponding to arrays and columns to coefficients to be estimated
f.model <- factor(targets, levels = levels(factor(targets)))

mod <- model.matrix(~0 + f.model)

colnames(mod) <- levels(factor(targets))


#####  Fit linear model for each gene given a series of arrays
fit <- lmFit(data_combat, design = mod)



#####  Create matrix of possible comparisons
comb <- combn(levels(factor(targets)), 2)

#####  Get number of possible comparisons using the following formula:
#
# n!/((n-r)!(r!))
#
# n = the number of classes to compare
# r = the number of elements for single comparison
#
################################################################################

#####   Record the number of groups and comparisons
targetsNo <- length(levels(factor(targets)))
combNo <- factorial(targetsNo)/(factorial(targetsNo-2)*(factorial(2))) # n!/((n-r)!(r!))

contrasts <- NULL
contrastNames <- NULL


#####  Create string with possible contrasts
for (i in 1:combNo) {
    
    contrasts <- c(contrasts, paste(paste(comb[1,i], comb[2,i], sep="vs"),paste(comb[1,i], comb[2,i], sep="-"), sep="="))
    contrastNames[i] <- paste(comb[1,i], comb[2,i], sep=" vs ")
}

contrasts <- paste(contrasts, collapse=", ")

#####  Create contrasts of interest
func = "makeContrasts"
arguments = paste(contrasts, "levels=mod",sep = ", ")

contrast.matrix <- eval(parse(text = paste(func, "(", arguments, ")", sep = "")))

#####  Fit contrasts to linear model
fitContrasts <- contrasts.fit(fit, contrast.matrix)

#####  Apply empirical Bayes statistics
eb <- eBayes(fitContrasts)



################################################################################
#
#    INTEGRATIVE ANALYSIS (3):
#
#    Results annotation
#
################################################################################


#####  Set preferred p-value and log2 fold-change threshold for calling differentially expressed genes
pThreshold = 0.0001
lfcThreshold = 2


#===============================================================================
#     Annotate probesets and write results into a file
#===============================================================================

#####  ... for each comparison
for (i in 1:ncol(eb$p.value) ) {
    
    topGenes <- topTable(eb, coef = colnames(eb)[i], adjust = "BH",  sort.by = "P", number = nrow(data_combat), genelist = rownames(data_combat))
    
    ##### Retrieve probesets annotation information
    probeid <- topGenes$ID
    
    SYMBOL <- as.character(unlist(lapply(mget(probeid,env=hgu133plus2SYMBOL), function (symbol) { return(paste(symbol,collapse="; ")) } )))
    NAME <- as.character(unlist(lapply(mget(probeid,env=hgu133plus2GENENAME), function (name) { return(paste(name,collapse="; ")) } )))
    CHR <- as.character(unlist(lapply(mget(probeid,env=hgu133plus2CHR), function (Chr) { return(paste(Chr,collapse="; ")) } )))
    MAP <- as.character(unlist(lapply(mget(probeid,env=hgu133plus2MAP), function (MAP) { return(paste(MAP,collapse="; ")) } )))
    
    ##### Merge probes annotation with statistics
    Annot <- data.frame(probeid, CHR, MAP, SYMBOL, NAME, row.names = NULL)
    Annot <- merge(Annot, topGenes, by.x = "probeid", by.y = "ID" )
    
    ##### Write annotated results into a file
    write.csv(Annot, file=paste("Integ_", colnames(eb$p.value)[i], "_topTable.csv", sep=""), row.names=FALSE)
}


#===============================================================================
#     ... only differnetially expressed genes
#===============================================================================

#####  Record differentially expressed genes
topGenes.index <- NULL

for (i in 1:ncol(eb$p.value) ) {
    
    topGenes <- topTable(eb, coef = colnames(eb)[i], adjust = "BH",  sort.by = "P", p.value = pThreshold, lfc = lfcThreshold, number = nrow(data_combat), genelist = rownames(data_combat))
    
    ##### Retrieve probesets annotation information
    probeid <- topGenes$ID
    
    SYMBOL <- as.character(unlist(lapply(mget(probeid,env=hgu133plus2SYMBOL), function (symbol) { return(paste(symbol,collapse="; ")) } )))
    NAME <- as.character(unlist(lapply(mget(probeid,env=hgu133plus2GENENAME), function (name) { return(paste(name,collapse="; ")) } )))
    CHR <- as.character(unlist(lapply(mget(probeid,env=hgu133plus2CHR), function (Chr) { return(paste(Chr,collapse="; ")) } )))
    MAP <- as.character(unlist(lapply(mget(probeid,env=hgu133plus2MAP), function (MAP) { return(paste(MAP,collapse="; ")) } )))
    
    ##### Merge probes annotation with statistics
    Annot <- data.frame(probeid, CHR, MAP, SYMBOL, NAME, row.names = NULL)
    Annot <- merge(Annot, topGenes, by.x = "probeid", by.y = "ID" )
    
    ##### Write annotated results into a file
    write.csv(Annot, file=paste("Integ_", colnames(eb$p.value)[i], "_DE.csv", sep=""), row.names=FALSE)
    
    #####  Record differentially expressed genes
    topGenes.index <- c(topGenes.index,rownames(topGenes))
}


################################################################################
#
#    INTEGRATIVE ANALYSIS (3):
#
#    Results visualisation
#
################################################################################


#===============================================================================
#     P-values histograms
#===============================================================================

for (i in 1:ncol(eb$p.value) ) {
    
    pdf(file = paste("Integ_", colnames(eb$p.value)[i], "_P_hist.pdf", sep = ""))
    histogram <- hist(eb$p.value[,i], breaks = seq(0,1, by = 0.01), main = contrastNames[i], xlab = "p-value")
    abline(v=0.05,col="red")
    dev.off()
}


#===============================================================================
#     Volcano plots
#===============================================================================

#####  Generate volcano plots of log2 fold-changes versus significance (adjusted p-values, label top 10 genes)
for (i in 1:ncol(eb$p.value) ) {
    
    topGenes <- topTable(eb, coef = colnames(eb)[i], adjust = "BH", sort.by = "none", number = nrow(data_combat))
    
    pdf(file = paste("Integ_", colnames(eb$p.value)[i], "_volcano_plot.pdf", sep = ""))
    plot(topGenes[,"logFC"], -log2(topGenes[,"adj.P.Val"]), pch = 16, cex = 0.5, xlab = "Log2 fold-change", ylab = "-log2(adjusted p-value)", main = contrastNames[i], col = "grey")
    
    #####  Highlight genes with logFC above specified threshold
    points(topGenes[abs(topGenes[, "logFC"]) > lfcThreshold, "logFC"], -log2(topGenes[abs(topGenes[, "logFC"]) > lfcThreshold, "adj.P.Val"]), cex = 0.5, pch = 16)
    abline(h = -log2(pThreshold), col = "red", lty = 2)
    
    #####  Label top 10 most significant genes
    ord <- order(-log2(topGenes[, "adj.P.Val"]), decreasing = TRUE)
    top10 <- ord[1:10]
    text(topGenes[top10, "logFC"], -log2(topGenes[top10, "adj.P.Val"]), labels = rownames(data_combat[top10,]), cex = 0.6, col = "blue")
    dev.off()
}


#===============================================================================
#     Hierarchical clustering
#===============================================================================

#####  Select top differentially expressed genes
topGenes <- data_combat[unique(topGenes.index),]


#####  Compute distance matrix and hierarchical clustering
hr <- hclust(as.dist(1-cor(t(topGenes), method = "pearson")), method = "ward.D2")
hc <- hclust(as.dist(dist(t(topGenes), method = "euclidean")), method = "ward.D2")


#####  Generate dendrogram
pdf("Integ_DE_dendrogram.pdf", width = 0.2*ncol(data)+2, height = 6, pointsize = 12)
par(mar = c(2,5,2,0))
plot(hc, xlab="", labels=paste(colnames(data_combat), targets, sep="       "), hang = -1, main="")
dev.off()


#####  ...heatmap (blue-red colour scale)
pdf("Integ_DE_heatmap_blue_red.pdf", width = 6, height = 10, pointsize = 12)
heatmap.2(as.matrix(topGenes), Rowv = as.dendrogram(hr), Colv=as.dendrogram(hc), col = colorRampPalette(c("blue","white","red"))(100), scale = "row", ColSideColors = targets.colour[[2]], margins = c(2, 6), labRow = "", labCol = "", trace = "none", key = TRUE)

#####  Add the legend
legend("topright", legend = levels(factor(targets)), fill = targets.colour[[1]], box.col = "transparent", cex=1.2)
dev.off()


#####  ...heatmap (green-red colour scale)
pdf("Integ_DE_heatmap_green_red.pdf", width = 6, height = 10, pointsize = 12)
heatmap.2(as.matrix(topGenes), Rowv = as.dendrogram(hr), Colv=as.dendrogram(hc), col = greenred(75), scale = "row", ColSideColors = targets.colour[[2]], margins = c(2, 6), labRow = "", labCol = "", trace = "none", key = TRUE)

#####  Add the legend
legend("topright", legend = levels(factor(targets)), fill = targets.colour[[1]], box.col = "transparent", cex=1.2)
dev.off()



##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

