# Functions
require(tidyverse)
require(Seurat)
require(colorblindr)
require(scrunchy)
# Function to make matrix and filter mtx

# cnt_path is path to the umi_tools output file
# bc_path is the path the the filtered barcode tsv file from 10x cellranger 3.0.0 output
# mtx_path is the path to the new haircut matrix directory
# filter_path is the path to the new filtered haircut matrix directory

make_mtx <- function(cnt_path, mtx_path, bc_path, filter_path) {
        #Convert UMI-tools output to matrix format
        umitools_to_mtx(count_file = cnt_path,
                        output_path = mtx_path)
        
        #Filter repair matrix for cell barcodes from cellranger 3.0.0 output
        filter_matrix(matrix_path = mtx_path,
                      barcodes_path = bc_path,
                      output_path = filter_path)
}

# Function to make and initialize a Seurat object

# mrna_path is the path to the filtered mrna matrix from cellranger 3 output
# haircut_path is the path to the filtered hiarcut matrix made through make_mtx above

make_seurat <- function(mrna_path, haircut_path){
        #read data
        mrna <- Read10X(data.dir = mrna_path)
        haircut <- Read10X(data.dir = haircut_path, gene.column = 1)
        
        #re-order repair matrix to match mrna matrix
        col.order <- colnames(mrna)
        haircut <- haircut[, col.order]
        
        #create seurat object and log normalize
        s <- CreateSeuratObject(counts = mrna)
        s[["repair"]] <- CreateAssayObject(counts = haircut)
        s = NormalizeData(object = s, assay = "RNA", method = "LogNormalize")
        s = NormalizeData(object = s, assay = "repair", method = "LogNormalize")
        
        # Calculating percent mitochondrial reads per cell 
        s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
        
        #Return seurat object
        s
}



# Function to scale data, classify PBMCs using TransferData from Seurat

# object is a seurat object and reference is a seurat object containing reference data
# and celltypes

analyze_seurat <- function(object, reference){
        
        # Scale data and run PCA
        object <- ScaleData(object, verbose = FALSE)
        object <- FindVariableFeatures(object, nfeatures = 2000)
        object <- RunPCA(object, npcs = 20, verbose = FALSE)
        # t-SNE and Clustering
        object <- RunTSNE(object, reduction = "pca", dims = 1:20)
        object <- RunUMAP(object, reduction = 'pca', dims = 1:20, spread = 2, min.dist = 1.5)
        object <- FindNeighbors(object, reduction = "pca", dims = 1:20)
        object <- FindClusters(object, resolution = 0.5)
        
        #Find anchors and transfer celltype labels from reference data
        marker.anchors <- FindTransferAnchors(reference = reference, query = object, 
                                              dims = 1:30)
        predictions <- TransferData(anchorset = marker.anchors, refdata = Idents(reference), 
                                    dims = 1:30)
        object <- AddMetaData(object, metadata = predictions)

        object
}

## Add simplified celltypes to object

# object is seurat object
# df is a dataframe of old names and new names
# var = variable to rename
add_simple_celltype <- function(object, df, var = "predicted.id"){
        
        ct <- rownames_to_column(FetchData(object, vars = var), "cell_id") %>% 
                left_join(df)
        
        object$celltype <- ct$celltype
        object
        
}