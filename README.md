# Single-Cell RNA-seq Analysis of NSCLC Tumor Cells

 **Dataset**: 20k dissociated tumor cells from NSCLC patients (10x Genomics)  
 **Tool used**: Seurat in R  
 **Goal**: Identify cellular heterogeneity and marker genes

## Contents
- `NSCLC_analysis.R`: Seurat pipeline
- `figures/`: All saved plots
- `output/`: Marker gene tables, Seurat object

## Key Steps
- Data import and quality control
- Normalization and variable feature detection
- PCA and UMAP-based clustering
- Marker gene identification

## Example Results
![UMAP](figures/umap_clusters.png)
![Heatmap](figures/marker_heatmap.png)

## Takeaways
- Identified 5 transcriptionally distinct clusters
- Top markers: ITGA1, TRGC2, LILRB4, FCGR2A
- Strong evidence of epithelial and immune cell populations

## Reference
[10x Genomics dataset](https://www.10xgenomics.com/datasets/20k-mixture-of-nsclc-dtcs-from-7-donors-3-v3-1-with-intronic-reads-3-1-standard)

