# Single-Cell RNA-seq Analysis of NSCLC Tumor Cells

 **Author:** Riya Dua
 **Dataset**: 20k dissociated tumor cells from NSCLC patients (10x Genomics)  
 **Tool used**: Seurat in R  
 **Goal**: Identify cellular heterogeneity and marker genes

## Objective:
To identify transcriptionally distinct subpopulations within NSCLC tumor samples using single-cell RNA sequencing (scRNA-seq), and characterize them based on their marker gene expression.


## Repository structure:
- `NSCLC_analysis.R`: Seurat pipeline
- `figures/`: All saved output plots
- `output/`: Marker gene tables, Seurat object (optional)

## Key Steps:
1.  Data import and quality control
2.  Normalization and variable feature detection
3.  PCA and UMAP-based clustering
4.  Marker gene identification

## Analysis Pipeline

| Step | Description |
|------|-------------|
| **1. Data Loading** | Read 10x `.h5` matrix using `Read10X_h5()` |
| **2. QC Filtering** | Removed low-quality cells (<200 or >2500 genes, >5% mitochondrial content) |
| **3. Normalization** | Log normalization of raw UMI counts |
| **4. Feature Selection** | Identified top 2,000 highly variable genes |
| **5. Scaling** | Centered and scaled data, regressing out unwanted variation |
| **6. PCA** | Performed linear dimensionality reduction |
| **7. Clustering** | Constructed SNN graph and identified clusters at multiple resolutions |
| **8. UMAP** | Visualized clusters in 2D space |
| **9. Marker Detection** | Identified top marker genes for each cluster |

## 📊 Key Visualizations

### 🔹 Quality Control
<p align="center"><img src="figures/QC_violinplot.png" width="500"/></p>

### 🔹 UMAP Clustering
<p align="center"><img src="figures/umap_clusters_res0.5.png" width="500"/></p>

### 🔹 Marker Gene Heatmap
<p align="center"><img src="figures/marker_heatmap.png" width="500"/></p>


## Takeaways
- Identified 5 transcriptionally distinct clusters
- Top markers: ITGA1, TRGC2, LILRB4, FCGR2A
- Strong evidence of epithelial and immune cell populations

## Referrences:
- [Seurat Documentation](https://satijalab.org/seurat/)
- [10x Genomics dataset](https://www.10xgenomics.com/datasets/20k-mixture-of-nsclc-dtcs-from-7-donors-3-v3-1-with-intronic-reads-3-1-standard)

