## PBMC make haircut matrices and filter matrices
## Make haircut matrices and filter matrices

library(Seurat)
library(tidyverse)
source("scripts/functions.R")

make_mtx("../data/pbmc/pbmc1/umitools_count.tsv.gz",
         "../data/pbmc/pbmc1/pbmc1_haircut",
         "../data/pbmc/pbmc1/filtered_feature_bc_matrix/barcodes.tsv.gz",
         "../data/pbmc/pbmc1/filtered_pbmc_haircut")
# there are 4009 barcodes remaining in the filtered data


make_mtx("../data/pbmc/pbmc2/umitools_counts.tsv.gz",
         "../data/pbmc/pbmc2/pbmc2_haircit",
         "../data/pbmc/pbmc2/filtered_feature_bc_matrix/barcodes.tsv.gz",
         "../data/pbmc/pbmc2/filtered_pbmc_haircut")
# there are 4023 barcodes remaining in the filtered data

make_mtx("../data/pbmc/pbmc3/umitools_counts.tsv.gz",
         "../data/pbmc/pbmc3/pbmc3_haircut",
         "../data/pbmc/pbmc3/filtered_feature_bc_matrix/barcodes.tsv.gz",
         "../data/pbmc/pbmc3/filtered_pbmc_haircut")
# there are 5140 barcodes remaining in the filtered data

# For dilution data

make_mtx("../data/pbmc/pbmc15/umitools_counts.tsv.gz",
         "../data/pbmc/pbmc15/pbmc15_haircut",
         "../data/pbmc/pbmc15/pbmc_15/barcodes.tsv.gz",
         "../data/pbmc/pbmc15/filtered_pbmc_haircut")
# there are 810 barcodes remaining in the filtered data

make_mtx("../data/pbmc/pbmc30/umitools_counts.tsv.gz",
         "../data/pbmc/pbmc30/pbmc30_haircut",
         "../data/pbmc/pbmc30/pbmc_30/barcodes.tsv.gz",
         "../data/pbmc/pbmc30/filtered_pbmc_haircut")
# there are 686 barcodes remaining in the filtered data

make_mtx("../data/pbmc/pbmc60/umitools_counts.tsv.gz",
         "../data/pbmc/pbmc60/pbmc60_haircut",
         "../data/pbmc/pbmc60/pbmc_60/barcodes.tsv.gz",
         "../data/pbmc/pbmc60/filtered_pbmc_haircut")
# there are 1383 barcodes remaining in the filtered data


## Making Seurat objects

pbmc1 <- make_seurat("../data/pbmc/pbmc1/filtered_feature_bc_matrix/",
                     "../data/pbmc/pbmc1/filtered_pbmc_haircut")

pbmc2 <- make_seurat("../data/pbmc/pbmc2/filtered_feature_bc_matrix/",
                     "../data/pbmc/pbmc2/filtered_pbmc_haircut/")

pbmc3 <- make_seurat("../data/pbmc/pbmc3/filtered_feature_bc_matrix/",
                     "../data/pbmc/pbmc3/filtered_pbmc_haircut/")


## For dilution/timecourse
pbmc15 <- make_seurat("../data/pbmc/pbmc15/pbmc_15/",
                     "../data/pbmc/pbmc15/filtered_pbmc_haircut/")

pbmc30 <- make_seurat("../data/pbmc/pbmc30/pbmc_30/",
                      "../data/pbmc/pbmc30/filtered_pbmc_haircut/")

pbmc60 <- make_seurat("../data/pbmc/pbmc60/pbmc_60/",
                      "../data/pbmc/pbmc60/filtered_pbmc_haircut/")
#### Determining cutoffs for number of features and percent mitochondrial reads per sample

#Make violin plot of features and cutoff outliers above and below feature cutoff
# and below % mitochondrial reads
VlnPlot(pbmc1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Decide on features and percent.mt cutoffs by looking at VlnPlot above
pbmc1 <- subset(pbmc1, subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 2000 & 
                        percent.mt < 15)

# Repeat for PBMC2 replicate
VlnPlot(pbmc2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Decide on features and percent.mt cutoffs by looking at VlnPlot above
pbmc2 <- subset(pbmc2, subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 2500 & 
                        percent.mt < 15)

# Repeat for PBMC3 replicate
VlnPlot(pbmc3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Decide on features and percent.mt cutoffs by looking at VlnPlot above
pbmc3 <- subset(pbmc3, subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 2000 & 
                        percent.mt < 15)

# Repeat for PBMC15 
VlnPlot(pbmc15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Decide on features and percent.mt cutoffs by looking at VlnPlot above
pbmc15 <- subset(pbmc15, subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 1000 & 
                        percent.mt < 15)

# Repeat for PBMC30 
VlnPlot(pbmc30, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Decide on features and percent.mt cutoffs by looking at VlnPlot above
pbmc30 <- subset(pbmc30, subset = nFeature_RNA > 200 & 
                         nFeature_RNA < 9000 & 
                         percent.mt < 15)

# Repeat for PBMC60 
VlnPlot(pbmc60, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Decide on features and percent.mt cutoffs by looking at VlnPlot above
pbmc60 <- subset(pbmc60, subset = nFeature_RNA > 200 & 
                         nFeature_RNA < 1500 & 
                         percent.mt < 15)

#load reference data
load("../data/other/pbmc_3k_reference.Rdata")

# Find variable features, calculate clusters, UMAP, TSNE
# Assign celltypes from reference single cell dataset

pbmc1 <- analyze_seurat(pbmc1, pbmc)
pbmc2 <- analyze_seurat(pbmc2, pbmc)
pbmc3 <- analyze_seurat(pbmc3, pbmc)
pbmc15 <- analyze_seurat(pbmc15, pbmc)
pbmc30 <- analyze_seurat(pbmc30, pbmc)
pbmc60 <- analyze_seurat(pbmc60, pbmc)

## Adding simplified celltype names to seurat objects
df <- tribble(~predicted.id, ~celltype,
              "CD8 T", "T",
              "CD14+ Mono", "Mono",
              "Naive CD4 T", "T",
              "NK", "NK",
              "Platelet", "Platelet",
              "B", "B",
              "DC",   "DC",        
              "FCGR3A+ Mono","Mono",
              "Memory CD4 T", "T")

pbmc1 <- add_simple_celltype(pbmc1, df)
pbmc2 <- add_simple_celltype(pbmc2, df)
pbmc3 <- add_simple_celltype(pbmc3, df)

pbmc15 <- add_simple_celltype(pbmc15, df)
pbmc30 <- add_simple_celltype(pbmc30, df)
pbmc60 <- add_simple_celltype(pbmc60, df)



# Save seurat objects
save(pbmc1, file = "../data/pbmc/seurat/pbmc1.seurat.Rdata")
save(pbmc2, file = "../data/pbmc/seurat/pbmc2.seurat.Rdata")
save(pbmc3, file = "../data/pbmc/seurat/pbmc3.seurat.Rdata")

save(pbmc15, file = "../data/pbmc/seurat/pbmc15.seurat.Rdata")
save(pbmc30, file = "../data/pbmc/seurat/pbmc30.seurat.Rdata")
save(pbmc60, file = "../data/pbmc/seurat/pbmc60.seurat.Rdata")

