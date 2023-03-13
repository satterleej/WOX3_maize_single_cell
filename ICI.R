library(psych)
library(dplyr)
library(Seurat)
library(data.table)
library(ggplot2)
library(MM)
library(RColorBrewer)
library(gplots)


#Index of cell identity

#Load meta data for computed ICI values, if ICI has already been calculated
#ICI.meta <- read.table(file = "/Users/jacksatterlee/Documents/Projects/narrowsheath/Ad_v_Ab_Timmermans/Timmermans_data_only_metadata.txt", row.names = 1, header = T)
#shoot.apex.ns.merged@meta.data <- ICI.meta

#Load in output of HTSeq-count
cpm <- read.table(file = "/Users/jacksatterlee/Documents/Projects/narrowsheath/Ad_v_Ab_Timmermans/cpm_Timmermans.txt", row.names = 1, header = T)

#Read in adaxial and abaxial identity genes
adaxial.genes <- row.names(read.table(file = "/Users/jacksatterlee/Documents/Projects/narrowsheath/Ad_v_Ab_Timmermans/adaxial_genes.txt", row.names = 1, header = T))
abaxial.genes <- row.names(read.table(file = "/Users/jacksatterlee/Documents/Projects/narrowsheath/Ad_v_Ab_Timmermans/abaxial_genes.txt", row.names = 1, header = T))

#Subset expression data by adaxial and abaxial identity genes with expression greater than 1 UMI. 
norm.exp.data <- apply(as.matrix(shoot.apex.ns.merged@assays$RNA@counts), 2, FUN = function(x){x/sum(x)})
row.names(norm.exp.data) <- gsub("gene:", "", row.names(norm.exp.data))
row.names(norm.exp.data) <- gsub("-F", "_F", row.names(norm.exp.data))

#Get number of cells
n_cell <- dim(norm.exp.data)[2]

#Get number of genes from adaxial and abaxial DE gene lists
n_gene_adaxial <- length(adaxial.genes)
n_gene_abaxial <- length(abaxial.genes)

#Get cpm data for DE genes from LCM-RNA-Seq
cpm.adaxial <- as.matrix(cpm[row.names(cpm) %in% adaxial.genes,])
cpm.abaxial <- as.matrix(cpm[row.names(cpm) %in% abaxial.genes,])

#Get normalized single cell expression data for DE genes from LCM-RNA-Seq
norm.exp.data.adaxial <- norm.exp.data[row.names(norm.exp.data) %in% adaxial.genes,]
norm.exp.data.abaxial <- norm.exp.data[row.names(norm.exp.data) %in% abaxial.genes,]

#Re-order genes in cpm and FC matrices in  alphanumeric order
cpm.adaxial <- cpm.adaxial[order(row.names(cpm.adaxial)),]
cpm.abaxial <- cpm.abaxial[order(row.names(cpm.abaxial)),]

#Re-order genes to match position in ordered cpm matrix
norm.exp.data.adaxial <- norm.exp.data.adaxial[order(row.names(norm.exp.data.adaxial)),]
norm.exp.data.abaxial <- norm.exp.data.abaxial[order(row.names(norm.exp.data.abaxial)),]

#Calculate spec scores from LCM RNA-Seq data
adaxial.spec.score <- vector(length = n_gene_adaxial)
for (i in 1:length(adaxial.genes)){
  adaxial.spec.score[i] <- ((mean(cpm.adaxial[i,1:2]) - mean(cpm.adaxial[i,3:4]))/(mean(cpm.adaxial[i,1:2]) + mean(cpm.adaxial[i,3:4])))
}

abaxial.spec.score <- vector(length = n_gene_abaxial)
for (i in 1:length(abaxial.genes)){
  abaxial.spec.score[i] <- ((mean(cpm.abaxial[i,3:4]) - mean(cpm.abaxial[i,1:2]))/(mean(cpm.abaxial[i,1:2]) + mean(cpm.abaxial[i,3:4])))
}

#Create matrix to store weighted gene expression data
adaxial <- matrix(nrow = n_gene_adaxial, ncol = n_cell, dimnames = list(row.names(adaxial.genes), colnames(norm.exp.data)))
abaxial <- matrix(nrow = n_gene_abaxial, ncol = n_cell, dimnames = list(row.names(abaxial.genes), colnames(norm.exp.data)))

#Create vectors to store the proportion of markers expressed in each cell
ratio.expressed.markers.adaxial <- vector(length = n_cell)
ratio.expressed.markers.abaxial <- vector(length = n_cell)
names(ratio.expressed.markers.adaxial) <- colnames(norm.exp.data)
names(ratio.expressed.markers.adaxial) <- colnames(norm.exp.data)

#Weight single cell gene expression based on spec scores
for (j in 1:n_cell){
  for (i in 1:n_gene_adaxial){
    adaxial[i,j] <- norm.exp.data.adaxial[i,j]*adaxial.spec.score[i]
  }
  ratio.expressed.markers.adaxial[j] <- length(which(norm.exp.data.adaxial[,j] != 0 & names(norm.exp.data.adaxial[,j]) %in% adaxial.genes))/n_gene_adaxial
}

for (j in 1:n_cell){
  for (i in 1:n_gene_abaxial){
    abaxial[i,j] <- norm.exp.data.abaxial[i,j]*abaxial.spec.score[i]
  }
  ratio.expressed.markers.abaxial[j] <- length(which(norm.exp.data.abaxial[,j] != 0 & names(norm.exp.data.abaxial[,j]) %in% abaxial.genes))/n_gene_abaxial
}

#Sum weighted expression data and normalize by number of total marker genes
adaxial.sum <- colSums(adaxial)/n_gene_adaxial
abaxial.sum <- colSums(abaxial)/n_gene_abaxial

#Further weight data by proportion of markers expressed in each cell
ICI.adaxial <- adaxial.sum*ratio.expressed.markers.adaxial
ICI.abaxial <- abaxial.sum*ratio.expressed.markers.abaxial

##############
#Random permutation analysis
##############

#Create matrix to store weighted gene expression data
adaxial <- matrix(nrow = n_gene_adaxial, ncol = n_cell, dimnames = list(row.names(adaxial.genes), colnames(norm.exp.data)))
abaxial <- matrix(nrow = n_gene_abaxial, ncol = n_cell, dimnames = list(row.names(abaxial.genes), colnames(norm.exp.data)))

n_iter <- 1000
ICI.adaxial.random <- matrix(nrow = n_cell, ncol = n_iter)

for (x in 1:n_iter){
  norm.exp.data.adaxial <- norm.exp.data[sample(1:dim(norm.exp.data)[1], n_gene_adaxial, replace=FALSE),] #Randomly select rows as adaxial marker genes
  for (j in 1:n_cell){
    for (i in 1:n_gene_adaxial){
      adaxial[i,j] <- norm.exp.data.adaxial[i,j]*adaxial.spec.score[i]
    }
    ratio.expressed.markers.adaxial[j] <- length(which(norm.exp.data.adaxial[,j] != 0))/n_gene_adaxial
  }
  ICI.adaxial.random[,x] <- (colSums(adaxial)/n_gene_adaxial)*ratio.expressed.markers.adaxial
  print(paste("Done with iteration",x,"for adaxial permutation"))
  
}

n_iter <- 1000
ICI.abaxial.random <- matrix(nrow = n_cell, ncol = n_iter)

for (x in 1:n_iter){
  norm.exp.data.abaxial <- norm.exp.data[sample(1:dim(norm.exp.data)[1], n_gene_abaxial, replace=FALSE),] #Randomly select rows as abaxial marker genes
  for (j in 1:n_cell){
    for (i in 1:n_gene_abaxial){
      abaxial[i,j] <- norm.exp.data.abaxial[i,j]*abaxial.spec.score[i]
    }
    ratio.expressed.markers.abaxial[j] <- length(which(norm.exp.data.abaxial[,j] != 0))/n_gene_abaxial
  }
  ICI.abaxial.random[,x] <- (colSums(abaxial)/n_gene_abaxial)*ratio.expressed.markers.abaxial
  print(paste("Done with iteration",x,"for abaxial permutation"))
  
}

####################
####################

#Normalize measured and random ICI scores such that they sum to 1
ICI.adaxial.norm <- ICI.adaxial/(ICI.adaxial + ICI.abaxial)
ICI.abaxial.norm <- ICI.abaxial/(ICI.abaxial + ICI.adaxial)
ICI.adaxial.random.norm <- ICI.adaxial.random/(ICI.adaxial.random + ICI.abaxial.random)
ICI.abaxial.random.norm <- ICI.abaxial.random/(ICI.abaxial.random + ICI.adaxial.random)


#Identify top 5% threshold of randomly generated ICI values for adaxial and abaxial ICIs
ICI.adaxial.random.norm.quantile <- vector(length = n_cell)
for (i in 1:n_cell){
  ICI.adaxial.random.norm.quantile[i] <- quantile(ICI.adaxial.random.norm[i,], probs = 0.05, na.rm = T)
}

ICI.abaxial.random.norm.quantile <- vector(length = n_cell)
for (i in 1:n_cell){
  ICI.abaxial.random.norm.quantile[i] <- quantile(ICI.abaxial.random.norm[i,], probs = 0.05, na.rm = T)
}

ICI.adaxial.identity <- ifelse(ICI.adaxial.norm >= ICI.adaxial.random.norm.quantile, "adaxial", "not adaxial")
ICI.abaxial.identity <- ifelse(ICI.abaxial.norm >= ICI.abaxial.random.norm.quantile, "abaxial", "not abaxial")

ICI.adaxial.identity.score <- ifelse(ICI.adaxial.identity == "adaxial", ICI.adaxial.norm, 0)
ICI.abaxial.identity.score <- ifelse(ICI.abaxial.identity == "abaxial", ICI.abaxial.norm, 0)

shoot.apex.ns.merged@meta.data$ICI.adaxial.identity.score <- ICI.adaxial.identity.score
shoot.apex.ns.merged@meta.data$ICI.abaxial.identity.score <- ICI.abaxial.identity.score

shoot.apex.ns.merged@meta.data$ICI.adaxial <- (ICI.adaxial.identity.score) - (ICI.abaxial.identity.score)
shoot.apex.ns.merged@meta.data$ICI.abaxial <- (ICI.abaxial.identity.score) - (ICI.adaxial.identity.score)

#Assign cells as either adaxial, abaxial, both, or neither
shoot.apex.ns.merged@meta.data$identity <- ifelse(shoot.apex.ns.merged@meta.data$ICI.adaxial.identity.score > 0 & shoot.apex.ns.merged@meta.data$ICI.abaxial.identity.score == 0, "adaxial",
                                                  ifelse(shoot.apex.ns.merged@meta.data$ICI.abaxial.identity.score > 0 & shoot.apex.ns.merged@meta.data$ICI.adaxial.identity.score == 0, "abaxial",
                                                         ifelse(shoot.apex.ns.merged@meta.data$ICI.adaxial.identity.score > 0 & shoot.apex.ns.merged@meta.data$ICI.abaxial.identity.score > 0, "both","neither")))

#Plot cells according to ad/ab identity
shoot.apex.ns.merged@meta.data$UMAP_1 <- shoot.apex.ns.merged@reductions$umap@cell.embeddings[,1]
shoot.apex.ns.merged@meta.data$UMAP_2 <- shoot.apex.ns.merged@reductions$umap@cell.embeddings[,2]

ggplot(data = shoot.apex.ns.merged@meta.data, aes(x = UMAP_1, y = UMAP_2, color = identity)) + geom_point(size = 0.1) +
  theme_classic() +
  scale_color_manual(values=c("goldenrod1", "cadetblue2", "grey"))
