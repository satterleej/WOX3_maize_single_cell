library(Seurat)
library(data.table)
library(ggplot2)
library(MM)
library(dplyr)

#Read 10X files output into R
shoot.apex.p6.ns.rep1 <- Read10X(data.dir = "~/Documents/Projects/narrowsheath/expression_matrices/NS_scRNAseq_outs_NS2/filtered_feature_bc_matrix_ns_scrnaseq_rep1", gene.column = 1)
shoot.apex.p6.WT.rep1 <- Read10X(data.dir = "~/Documents/Projects/narrowsheath/expression_matrices/NS_scRNAseq_outs_NS2/filtered_feature_bc_matrix_WT_scrnaseq_rep1", gene.column = 1)
shoot.apex.p6.ns.rep2 <- Read10X(data.dir = "~/Documents/Projects/narrowsheath/expression_matrices/NS_scRNAseq_outs_NS2/filtered_feature_bc_matrix_ns_scrnaseq_rep2", gene.column = 1)
shoot.apex.p6.WT.rep2 <- Read10X(data.dir = "~/Documents/Projects/narrowsheath/expression_matrices/NS_scRNAseq_outs_NS2/filtered_feature_bc_matrix_WT_scrnaseq_rep2", gene.column = 1)

#Create Seurat objects for each sample
shoot.apex.p6.ns.rep1 <- CreateSeuratObject(counts = shoot.apex.p6.ns.rep1, project = "shoot.apex.ns.rep1", min.cells = 0, min.features = 0)
shoot.apex.p6.WT.rep1 <- CreateSeuratObject(counts = shoot.apex.p6.WT.rep1, project = "shoot.apex.wt.rep1", min.cells = 0, min.features = 0)
shoot.apex.p6.ns.rep2 <- CreateSeuratObject(counts = shoot.apex.p6.ns.rep2, project = "shoot.apex.ns.rep2", min.cells = 0, min.features = 0)
shoot.apex.p6.WT.rep2 <- CreateSeuratObject(counts = shoot.apex.p6.WT.rep2, project = "shoot.apex.wt.rep2", min.cells = 0, min.features = 0)

#Load in a list of mitochondrial genes
mito.genes <- read.csv("~/Documents/Projects/narrowsheath/mitochondrial_genes.txt", row.names = 1)

#Calculate the percentage of mitochondrial genes expressed in each cell
shoot.apex.p6.ns.rep1[["percent.mt"]] <- PercentageFeatureSet(shoot.apex.p6.ns.rep1, features = shoot.apex.p6.ns.rep1@assays$RNA@counts@Dimnames[[1]][as.factor(row.names(mito.genes))])
shoot.apex.p6.WT.rep1[["percent.mt"]] <- PercentageFeatureSet(shoot.apex.p6.WT.rep1, features = shoot.apex.p6.WT.rep1@assays$RNA@counts@Dimnames[[1]][as.factor(row.names(mito.genes))])
shoot.apex.p6.ns.rep2[["percent.mt"]] <- PercentageFeatureSet(shoot.apex.p6.ns.rep2, features = shoot.apex.p6.ns.rep2@assays$RNA@counts@Dimnames[[1]][as.factor(row.names(mito.genes))])
shoot.apex.p6.WT.rep2[["percent.mt"]] <- PercentageFeatureSet(shoot.apex.p6.WT.rep2, features = shoot.apex.p6.WT.rep2@assays$RNA@counts@Dimnames[[1]][as.factor(row.names(mito.genes))])

#Filter cells based on number of genes expressed, number of unique transcripts, and percent mitochondrial genes expressed
VlnPlot(shoot.apex.p6.ns.rep1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
VlnPlot(shoot.apex.p6.WT.rep1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
VlnPlot(shoot.apex.p6.ns.rep2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
VlnPlot(shoot.apex.p6.WT.rep2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

shoot.apex.p6.ns.rep1 <- subset(shoot.apex.p6.ns.rep1, subset = percent.mt < 1 & nFeature_RNA > 3000 & nCount_RNA < 1e+05)
shoot.apex.p6.WT.rep1 <- subset(shoot.apex.p6.WT.rep1, subset = percent.mt < 1 & nFeature_RNA > 3000 & nCount_RNA < 1e+05)
shoot.apex.p6.ns.rep2 <- subset(shoot.apex.p6.ns.rep2, subset = percent.mt < 1 & nFeature_RNA > 3000 & nCount_RNA < 1e+05)
shoot.apex.p6.WT.rep2 <- subset(shoot.apex.p6.WT.rep2, subset = percent.mt < 1 & nFeature_RNA > 3000 & nCount_RNA < 1e+05)

#Merge the datasets from individual samples
shoot.apex.ns.merged <- merge(x = shoot.apex.p6.ns.rep1, y = c(shoot.apex.p6.WT.rep1, shoot.apex.p6.ns.rep2, shoot.apex.p6.WT.rep2),
                              add.cell.ids = c("ns", "WT", "ns", "WT"),
                              project = "shoot.apex.ns")

#Run data transformation, PCA, and k-nearest neighbors hierarchical clustering
shoot.apex.ns.merged <- SCTransform(shoot.apex.ns.merged)
shoot.apex.ns.merged <- RunPCA(shoot.apex.ns.merged, verbose = FALSE)
shoot.apex.ns.merged <- FindNeighbors(shoot.apex.ns.merged, dims = 1:25)
shoot.apex.ns.merged <- FindClusters(shoot.apex.ns.merged, resolution = 1)

#Identify cells expressing at least one WOX3 gene (WOX3a/b, NS1)
wox3.expression.summary <- shoot.apex.ns.merged@assays$RNA@counts[c("gene:GRMZM2G069028","gene:GRMZM2G140083","gene:GRMZM2G122537","gene:Zm00001d052598"),]
wox3.expression.summary <- as.array(wox3.expression.summary)
wox3.expression.summary.colsums <- colSums(wox3.expression.summary)
wox3.expressing.cells <- names(wox3.expression.summary.colsums[wox3.expression.summary.colsums > 0])

#Add categorical variable to meta.data for wox3 expression
shoot.apex.ns.merged@meta.data$expression <- ifelse(row.names(shoot.apex.ns.merged@meta.data) %in% wox3.expressing.cells, "positive", "negative")
#Add categorical variable to meta.data for genotype
shoot.apex.ns.merged@meta.data$genotype <- ifelse(shoot.apex.ns.merged@meta.data$orig.ident == "shoot.apex.wt.rep1" | shoot.apex.ns.merged@meta.data$orig.ident == "shoot.apex.wt.rep2" , "WT", "ns")
#Add categorical variable to meta.data for replicate
shoot.apex.ns.merged@meta.data$replicate <- ifelse(shoot.apex.ns.merged@meta.data$orig.ident == "shoot.apex.wt.rep1" | shoot.apex.ns.merged@meta.data$orig.ident == "shoot.apex.ns.rep1", "rep1", "rep2")


#Run UMAP for ordinations
shoot.apex.ns.merged <- RunUMAP(shoot.apex.ns.merged, dims = 1:25, n.neighbors = 25, min.dist = 0.01, spread = 1)
