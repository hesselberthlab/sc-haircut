# Functions
require(tidyverse)
require(Seurat)
require(colorblindr)
require(scrunchy)
require(ggplot2)
require(cowplot)

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
        s <- FindVariableFeatures(s, selection.method = "vst", nfeatures = 5000)
        
        #Return seurat object
        s
}



# Function to scale data, classify PBMCs using TransferData from Seurat

# object is a seurat object and reference is a seurat object containing reference data
# and celltypes

analyze_seurat <- function(object, reference){
        
        object <- NormalizeData(object, verbose = FALSE)
        object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 5000)
        
        #Find anchors and transfer celltype labels from reference data
        marker.anchors <- FindTransferAnchors(reference = reference, query = object, 
                                              dims = 1:30)
        predictions <- TransferData(anchorset = marker.anchors, refdata = Idents(reference), 
                                    dims = 1:30)
        object <- AddMetaData(object, metadata = predictions)
        
        # Scale data and run PCA
        object <- ScaleData(object, verbose = FALSE)
        object <- RunPCA(object, npcs = 20, verbose = FALSE)
        # t-SNE and Clustering
        object <- RunUMAP(object, reduction = "pca", dims = 1:20)
        object <- FindNeighbors(object, reduction = "pca", dims = 1:20)
        object <- FindClusters(object, resolution = 0.5)
        
        object <- RunUMAP(object, reduction = 'pca', dims = 1:10, spread = 2, min.dist = 1.5)
        
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


## Color palette
# colors
colors <- c(
        palette_OkabeIto_black,
        scales::brewer_pal(palette = "Paired")(12),
        scales::brewer_pal(palette = "Set1")(9),
        scales::brewer_pal(palette = "Set2")(8),
        scales::brewer_pal(palette = "Dark2")(8)
)
loupe_palette <- rev(scales::brewer_pal(palette = "RdGy")(11)[c(1:5, 7)])


# Bulk hairpin plot code
haircut_plot <-  function(data, x, y, 
                          col = "sample", 
                          x_lab = "Hairpin Position", 
                          y_lab = "Normalized Counts (count / total reads)", 
                          pal =  colors, 
                          xlim = c(0, 60),
                          point = F){
        ggplot(data, aes_string(x = x, y = y, color = col)) +
                geom_line(alpha=.7, size = .8) + 
                #        facet_wrap(~hairpin) + 
                geom_vline(xintercept = 44, linetype="dotted", color = "black", size = 1) +
                #geom_vline(xintercept = 10, linetype="dotted", color = "grey", size = 1) +
                theme_cowplot() + 
                theme(legend.position="top",
                      legend.title = element_blank()) +
                scale_color_manual(values = pal) +
                #        scale_colour_manual(values = rev(cols[1:2])) + 
                ylab(y_lab) + 
                xlab(x_lab) +
                xlim(xlim) + 
                if(point) geom_point()
}

## get hairpin coverage from seurat object

get_hairpin_coverage <- function(s, var = "celltype") {
        
        df <- FetchData(s, vars = var)
        
        colnames(df) <- str_replace(colnames(df), pattern = "-", replacement = "_")
        colnames(df) <- str_replace(colnames(df), pattern = "repair_", replacement = "")
        
        cnts = GetAssayData(object = s, slot = "counts", assay = 'repair')
        
        category <- df[, var]
        cell_ids <- rownames(df)
        cell_ids <- split(rownames(df), category)
        
        res <- purrr::map_dfr(
                cell_ids,
                function(ids) {
                        hairpin_info <- as.data.frame(rownames(cnts))
                        hairpin_info$count <- Matrix::rowSums(cnts[, ids, drop = FALSE])
                        hairpin_info
                },
                .id = "cell_id"
        )
        
        colnames(res) <- c("celltype", "hairpin_pos", "count")
        res <- separate(res, hairpin_pos, into = c("hairpin", "position")) %>%
                mutate(position = as.double(position))
        res <- as_tibble(res)
        
        df %>% group_by(celltype) %>% 
                tally() -> t
        full_join(res, t) %>%
                mutate(avg_count = count/n) -> res
        res
}


activity_plot <- function(df, vline = c(2, 4, 6), lab = NULL){
        lab <- enquo(lab)
        if(!is.null(lab)) { 
                p <- ggplot(df, aes(y = celltype, x = activity, color = celltype)) + 
                        ggbeeswarm::geom_quasirandom(aes(alpha = !!lab), 
                                                     size = 0.5, 
                                                     groupOnX = F,
                                                     dodge.width = 1) +
                        scale_color_manual(values = colors) + 
                        theme_cowplot() + 
                        theme(legend.position = 'none', 
                              strip.placement = "outside",
                              strip.background = element_blank(),
                        ) +
                        ylab(NULL) + 
                        geom_vline(xintercept = vline, color =  "#999999",
                                   alpha = 0.5, linetype = "dotted")
        }
        else {
        p <- ggplot(df, aes(y = celltype, x = activity, color = celltype)) + 
                ggbeeswarm::geom_quasirandom(size = 0.5, groupOnX = F) +
                scale_color_manual(values = colors) + 
                theme_cowplot() + 
                theme(legend.position = 'none', 
                      strip.placement = "outside",
                      strip.background = element_blank(),
                ) +
                ylab(NULL) + 
                geom_vline(xintercept = vline, color =  "#999999",
                           alpha = 0.5, linetype = "dotted") 
                }
        p
}

get_single_cell_df <- function(s, feat){
        df <- rownames_to_column(FetchData(pbmc1, feat), 
                                 "cell_id") %>% 
                mutate(cell_id = paste(cell_id, 1, sep = "_"))
        colnames(df) <- str_replace(colnames(df), "-", "_")
        colnames(df) <- str_replace(colnames(df), "repair_", "")
        df
}





