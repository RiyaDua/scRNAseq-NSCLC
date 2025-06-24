# Title: single cell analysis using seurat
# Author: Riya Dua
# Goal: Identify cellular heterogeneity and marker genes in Non-small cell lung cancer
# Dataset: 20k Mixture of Non-small cell lung cancer  (NSCLC) dissociated tumor cells (DTCs) from 7 donors, 3' v3.1 (without intronic reads)
#' Gene Expression libraries were generated from ~33,000 cells 
#' link: https://www.10xgenomics.com/datasets/20k-mixture-of-nsclc-dtcs-from-7-donors-3-v3-1-with-intronic-reads-3-1-standard


# 1. installing libraries ====
#install.packages("Seurat")
library("Seurat")
library(tidyverse)
library(hdf5r)

# 2. loading the dataset downaloaded from 10x Genomics ====
data<-Read10X_h5("Desktop/Independent projects/Seurat/data/20k_NSCLC_DTC_3p_nextgem_donor_1_count_sample_feature_bc_matrix.h5")
str(data)

counts<-data$`Gene Expression`
str(counts)


# 3. initializing seurat object ====
#' In our seurat object we provide counts (raw), name of the project
#' Minimum number of cells that we want -in which our gene (feature) is expressed
#' Minimum number of features (genes) that we want in all the cells -
#' so this keeps all those cells that have at least 200 features being expressed 

seurat.obj<-CreateSeuratObject(counts = counts, project = "NSCLC", min.cells = 3, min.features = 200)

# To checkout the slots created:
str(seurat.obj)
seurat.obj

metadata<-seurat.obj@meta.data


# 4. Quality Control ====
#' We have raw counts matrix and thus it needs to be filtered our of bad quality cells
#' we will look at total genes, cells, molecules and so on
#' Poor quality cells will have low number of genes and molecules
#' we will also look at % of mitochondrial genes - in low quality cells we see higher mitochondrial gene contamination

## 4.1. Mitochondrial percentage====
#' creating a column called percent.mt in our seurat object that will have calculated percentage of mitochondrial content
#' in each cell
seurat.obj[["percent.mt"]]<-PercentageFeatureSet(seurat.obj, pattern = "^MT-")


# visualizing these features in a violin plot
vln<-VlnPlot(object = seurat.obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
fs<-FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
ggsave("Desktop/Independent projects/Seurat/figures/QC_violinplot.png", plot = vln, width = 9, height = 4, dpi = 300)
ggsave("Desktop/Independent projects/Seurat/figures/Feature_scatter.png", plot = fs, dpi = 300)# width = 9, height = 4, dpi = 300)


## 4.2. Filtering out low quality cells ====
seurat.obj<-subset(seurat.obj, subset = nFeature_RNA>200 & nFeature_RNA <2500 & percent.mt <5)

## 4.3. Normalize data ====
seurat.obj<-NormalizeData(seurat.obj)


## 4.4. Identifying highly variable features ====
seurat.obj<-FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)


# top 10 most variable genes
top10<-head(VariableFeatures(seurat.obj), 10)

# plotting those features
plot1<-VariableFeaturePlot(seurat.obj)
plot_1<-LabelPoints(plot = plot1, points = top10, repel=TRUE)
ggsave("Desktop/Independent projects/Seurat/figures/Variable_features.png", plot = plot_1, dpi = 300)


## 4.5. scaling data ====
#' accounting for batch effect and other technical variations
genes<-rownames(seurat.obj)
seurat.obj<-ScaleData(seurat.obj, features = genes)
str(seurat.obj)


## 5. Linear dimensionality reduction ====
seurat.obj<-RunPCA(seurat.obj, features = VariableFeatures(seurat.obj))

# visualizing PCA results
print(seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)


# Loop through dims 1 to 5 and save each heatmap
for (i in 1:5) {
  file_name <- paste0("Desktop/Independent projects/Seurat/figures/pca_dimheatmap_PC", i, ".png")
  png(file_name, width = 1000, height = 800)
  DimHeatmap(seurat.obj, dims = i, cells = 500, balanced = TRUE)
  dev.off()
}


# determine dimensionality of the data
elb<-ElbowPlot(seurat.obj)
ggsave("Desktop/Independent projects/Seurat/figures/elbow_plot.png", plot = elb, dpi = 300)


# 6. Clustering ====
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:15)

# understanding resolution
seurat.obj <- FindClusters(seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(seurat.obj@meta.data)

# Based on UMAP and clustering, we observe ~6 distinct transcriptional clusters 
# which may correspond to different cell types/states in the NSCLC tumor
set.seed(123)
png("Desktop/Independent projects/Seurat/figures/cluster_dimplot.png", width = 1000, height = 800)
DimPlot(seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)
dev.off()

# setting identity of clusters
Idents(seurat.obj)
Idents(seurat.obj) <- "RNA_snn_res.0.1"
Idents(seurat.obj)

# 7. non-linear dimensionality reduction =====
seurat.obj <- RunUMAP(seurat.obj, dims = 1:15)

# 8. individual clusters ====
set.seed(123)
png("Desktop/Independent projects/Seurat/figures/ind_cluster_dimplot.png", width = 1000, height = 800)
DimPlot(seurat.obj, reduction = "umap")
dev.off()

## 8.1. Identifying markers for each cluster ====
#' This means we are running differential gene expression analysis between clusters to find marker genes
#' FindAllMarkers() runs DE for every cluster vs. every other cell one at a time
#' only.pos = TRUE returns only upregulated genes (basically genes that have positive log FC)
#' min.pct = 0.25 sets a threshold of significance that a gene must be expressed in 25% of cells in either cluster

markers <- FindAllMarkers(seurat.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)
write.csv(markers, "~/Desktop/Independent projects/Seurat/output/NSCLC_cluster_markers.csv", row.names = FALSE)

## 8.2. plotting top markers ====
top_marker<-markers %>%
  group_by(cluster)%>%
  top_n(n=2, wt = avg_log2FC) # selecting top 2 markers per cluster


#' In DoHeatmap, each row is a feature (gene) and each column represents a cell (not averaged)
#' Used to confirm if a gene is uniquely expressed in a particular cluster
#' Helps visualize how specific genes vary across all cells, especially within and between clusters 
png("Desktop/Independent projects/Seurat/figures/top_marker_genes.png", width = 1000, height = 800)
DoHeatmap(seurat.obj, features = top_marker$gene)
dev.off()


## 8.3. Feature plots (expression on UMAP) ====
#' to see expression of the top marker genes 
png("Desktop/Independent projects/Seurat/figures/feature_plot.png", width = 1000, height = 800)
FeaturePlot(seurat.obj, features = c("ITGA1", "TRGC2", "LILRB4", "FCGR2A")) 
dev.off()

## 8.4. DotPlot of marker expression across clusters ====
png("Desktop/Independent projects/Seurat/figures/marker_expression_per_cluster.png", width = 1000, height = 800)
DotPlot(seurat.obj, features = c("EPCAM", "VIM", "MKI67", "CDH1")) + RotatedAxis()
dev.off()

## 8.5. saving final seurat object =====
#saveRDS(seur.obj, "Desktop/Independent projects/Seurat/output/Seurat_object.rds")

# 9. Session Info ====
sessionInfo()
