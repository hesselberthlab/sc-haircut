## Make haircut matrices and filter matrices

library(scrunchy)
library(tidyverse)
library(Seurat)
source("scripts/functions.R")
#Convert UMI-tools output to matrix format

make_mtx(cnt_path = "../data/barnyard/umitools_counts.tsv.gz",
         mtx_path = "../data/barnyard/haircut_mtx",
         bc_path = "../data/barnyard/barnyard_cellranger3_out/barcodes.tsv.gz",
         filter_path = "../data/barnyard/filtered_haircut")

# there are 2299 barcodes remaining in the filtered data

## Making Seurat object
barnyard_seurat <- make_seurat(mrna_path = "../data/barnyard/barnyard_cellranger3_out",
                               haircut_path = "../data/barnyard/filtered_haircut")

## Find variable featuers, scale data, run PCA and TSNE and UMAP
barnyard_seurat <- FindVariableFeatures(object = barnyard_seurat, 
                                        selection.method = 'mean.var.plot', 
                                        mean.cutoff = c(0.0125, 3), 
                                        dispersion.cutoff = c(0.5, Inf))
barnyard_seurat <- ScaleData(object = barnyard_seurat, 
                             features = rownames(x = barnyard_seurat))

barnyard_seurat <- RunPCA(object = barnyard_seurat, 
                          features = VariableFeatures(object = barnyard_seurat), 
                          verbose = FALSE)

#use elbow plot to determine number of dims to bring forward
# ElbowPlot(object = barnyard_seurat) 

#Cluster and run UMAP 
barnyard_seurat <- FindNeighbors(object = barnyard_seurat, dims = 1:6)
barnyard_seurat <- FindClusters(object = barnyard_seurat, resolution = 0.6)
barnyard_seurat <- RunTSNE(object = barnyard_seurat, dims = 1:6)
barnyard_seurat <- RunUMAP(object = barnyard_seurat, dims = 1:6)

#Save seurat object for downstream use
save(barnyard_seurat, file = "../data/barnyard/barnyard.seurat.object.Rdata")


