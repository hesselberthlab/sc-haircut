# Supplementary figure 10
# Classifying cells with different reference data

source("scripts/functions.R")

## Classification of cells using mRNA alone, mRNA with repair, or repair alone

library(tidyverse)
library(Seurat)
source("scripts/functions.R")
library(cowplot)

load("../data/pbmc/seurat/pbmc1.seurat.Rdata")
load("../data/pbmc/seurat/pbmc2.seurat.Rdata")

## Transfer anchors using mRNA alone
Idents(pbmc1) <- pbmc1$celltype
marker.anchors <- FindTransferAnchors(reference = pbmc1, query = pbmc2, 
                                      dims = 1:30)

predictions <- TransferData(anchorset = marker.anchors, refdata = Idents(pbmc1), 
                            dims = 1:30)

pbmc2_new <- AddMetaData(pbmc2, metadata = predictions)

# Compare plots
plot_grid(DimPlot(pbmc2, reduction = 'umap', group.by = "celltype",
                  cols = colors) + 
                  ggtitle("Original classification"),
          DimPlot(pbmc2_new, reduction = 'umap', group.by = 'predicted.id',
                  cols = colors) + 
                  ggtitle("New classification mRNA only"))

# Get % of cells
orig <- rownames_to_column(FetchData(pbmc2, vars = "celltype"))
new <- rownames_to_column(FetchData(pbmc2_new, vars = "predicted.id"))

# rows are in the same order
sum(orig$rowname == new$rowname) / length(orig$rowname)

# 93% of cells are classified the same way
sum(orig$celltype == new$predicted.id) / length(orig$rowname)

colnames(orig) <- c("cell_id", "original")
colnames(new) <- c("cell_id", "new")

df <- full_join(orig, new)
df %>% mutate(same = ifelse(original == new, "same", 'different')) %>%
        group_by(original, same) %>%
        summarise(total = n()) %>%
        mutate(percent = total / sum(total)) %>%
        filter(same == "same") %>%
        ggplot(aes(x = original,y =  percent, fill = original)) +
        geom_bar(stat = 'identity') + 
        scale_fill_manual(values = colors)

## Make repair only object
pbmc1_repair <- CreateSeuratObject( counts =GetAssayData(pbmc1, 
                                                         assay = "repair", 
                                                         slot = "counts"))

pbmc1_repair <- NormalizeData(pbmc1_repair, verbose = FALSE)
pbmc1_repair <- FindVariableFeatures(pbmc1_repair, selection.method = "vst", nfeatures = 200)
pbmc1_repair <- ScaleData(pbmc1_repair, verbose = FALSE)
pbmc1_repair <- RunPCA(pbmc1_repair, npcs = 20, verbose = FALSE)
pbmc1_repair <- RunUMAP(pbmc1_repair, reduction = "pca", dims = 1:20)
pbmc1_repair <- FindNeighbors(pbmc1_repair, reduction = "pca", dims = 1:20)
pbmc1_repair <- FindClusters(pbmc1_repair, resolution = 0.5)

DimPlot(pbmc1_repair, reduction = 'umap', cols = colors)
FeaturePlot(pbmc1_repair, features = c("Uracil-45", "riboG-44", 
                                       "GU-45", "Abasic-45", "Abasic-46"),
            reduction = 'umap')

pbmc1_repair$predicted.id <- pbmc1$celltype
Idents(pbmc1_repair) <- pbmc1_repair$predicted.id
DimPlot(pbmc1_repair, reduction = 'umap', cols = colors)

# Filter PBMC2 to have only the repair substrates in PBMC1
rep_counts <- GetAssayData(pbmc2, assay = "repair", slot = "counts")
rep_counts <- rep_counts[rownames(pbmc1_repair), ]

pbmc2_repair <- CreateSeuratObject( counts = rep_counts)

pbmc2_repair <- NormalizeData(pbmc2_repair, verbose = FALSE)
pbmc2_repair <- FindVariableFeatures(pbmc2_repair, selection.method = "vst", nfeatures = 200)
pbmc2_repair <- ScaleData(pbmc2_repair, verbose = FALSE)
pbmc2_repair <- RunPCA(pbmc2_repair, npcs = 20, verbose = FALSE)
pbmc2_repair <- RunUMAP(pbmc2_repair, reduction = "pca", dims = 1:20)
pbmc2_repair <- FindNeighbors(pbmc2_repair, reduction = "pca", dims = 1:20)
pbmc2_repair <- FindClusters(pbmc2_repair, resolution = 0.5)

DimPlot(pbmc2_repair, reduction = 'umap', cols = colors)
pbmc2_repair$predicted.id <- pbmc2$celltype
Idents(pbmc2_repair) <- pbmc2_repair$predicted.id
DimPlot(pbmc2_repair, reduction = 'umap', cols = colors)

## Find anchors and transfer labels
marker.anchors <- FindTransferAnchors(reference = pbmc1_repair, query = pbmc2_repair, 
                                      dims = 1:30)

predictions <- TransferData(anchorset = marker.anchors, refdata = Idents(pbmc1_repair), 
                            dims = 1:30)

pbmc2_repair_new <- AddMetaData(pbmc2_repair, metadata = predictions)

# Compare plots
plot_grid(DimPlot(pbmc2_repair, reduction = 'umap', group.by = "predicted.id",
                  cols = colors) + 
                  ggtitle("Original classification"),
          DimPlot(pbmc2_repair_new, reduction = 'umap', group.by = 'predicted.id',
                  cols = colors) + 
                  ggtitle("New classification repair only"))

# Get % of cells
orig_rep <- rownames_to_column(FetchData(pbmc2_repair, vars = "predicted.id"))
new_rep <- rownames_to_column(FetchData(pbmc2_repair_new, vars = "predicted.id"))

# rows are in the same order
sum(orig_rep$rowname == new_rep$rowname) / length(orig_rep$rowname)

# 49%% of cells are classified the same way
sum(orig_rep$predicted.id == new_rep$predicted.id) / length(orig_rep$rowname)

colnames(orig_rep) <- c("cell_id", "original")
colnames(new_rep) <- c("cell_id", "new")

df <- full_join(orig_rep, new_rep)
df %>% mutate(same = ifelse(original == new, "same", 'different')) %>%
        group_by(original, same) %>%
        summarise(total = n()) %>%
        mutate(percent = total / sum(total)) %>%
        filter(same == "same") %>%
        ggplot(aes(x = original,y =  percent, fill = original)) +
        geom_bar(stat = 'identity') + 
        scale_fill_manual(values = colors)

## Combining repair and mRNA into one matrix

mrna <- GetAssayData(pbmc1, assay = "RNA", slot = 'counts')
repair <- GetAssayData(pbmc1, assay = 'repair', slot = 'counts')

pbmc1_comb <- CreateSeuratObject( counts = rbind(mrna, repair))

pbmc1_comb <- NormalizeData(pbmc1_comb, verbose = FALSE)
pbmc1_comb <- FindVariableFeatures(pbmc1_comb, selection.method = "vst", nfeatures = 5000)
pbmc1_comb <- ScaleData(pbmc1_comb, verbose = FALSE)
pbmc1_comb <- RunPCA(pbmc1_comb, npcs = 20, verbose = FALSE)
pbmc1_comb <- RunUMAP(pbmc1_comb, reduction = "pca", dims = 1:20)
pbmc1_comb <- FindNeighbors(pbmc1_comb, reduction = "pca", dims = 1:20)
pbmc1_comb <- FindClusters(pbmc1_comb, resolution = 0.5)

pbmc1_comb$predicted.id <- pbmc1$celltype
Idents(pbmc1_comb) <- pbmc1_comb$predicted.id
DimPlot(pbmc1_comb, reduction = 'umap', cols = colors)

# Uracil-45, Abasic-45, Abasic-46 all show up in the top 100 variable features
VariableFeatures(pbmc1_comb)

mrna <- GetAssayData(pbmc2, assay = "RNA", slot = 'counts')
repair <- GetAssayData(pbmc2, assay = 'repair', slot = 'counts')
repair <- repair[rownames(pbmc1_repair), ]

pbmc2_comb <- CreateSeuratObject( counts = rbind(mrna, repair))

pbmc2_comb <- NormalizeData(pbmc2_comb, verbose = FALSE)
pbmc2_comb <- FindVariableFeatures(pbmc2_comb, selection.method = "vst", nfeatures = 5000)
pbmc2_comb <- ScaleData(pbmc2_comb, verbose = FALSE)
pbmc2_comb <- RunPCA(pbmc2_comb, npcs = 20, verbose = FALSE)
pbmc2_comb <- RunUMAP(pbmc2_comb, reduction = "pca", dims = 1:20)
pbmc2_comb <- FindNeighbors(pbmc2_comb, reduction = "pca", dims = 1:20)
pbmc2_comb <- FindClusters(pbmc2_comb, resolution = 0.5)

pbmc2_comb$predicted.id <- pbmc2$celltype
Idents(pbmc2_comb) <- pbmc2_comb$predicted.id
DimPlot(pbmc2_comb, reduction = 'umap', cols = colors)

marker.anchors <- FindTransferAnchors(reference = pbmc1_comb, query = pbmc2_comb, 
                                      dims = 1:30)

predictions <- TransferData(anchorset = marker.anchors, refdata = Idents(pbmc1_comb), 
                            dims = 1:30)

pbmc2_comb_new <- AddMetaData(pbmc2_comb, metadata = predictions)

# Compare plots
plot_grid(DimPlot(pbmc2_comb, reduction = 'umap', group.by = "predicted.id",
                  cols = colors) + 
                  ggtitle("Original classification"),
          DimPlot(pbmc2_comb_new, reduction = 'umap', group.by = 'predicted.id',
                  cols = colors) + 
                  ggtitle("New classification mRNA and repair"))

# Get % of cells

orig_rep <- rownames_to_column(FetchData(pbmc2_comb, vars = "predicted.id"))
new_rep <- rownames_to_column(FetchData(pbmc2_comb_new, vars = "predicted.id"))

# rows are in the same order
sum(orig_rep$rowname == new_rep$rowname) / length(orig_rep$rowname)

# 93%% of cells are classified the same way
sum(orig_rep$predicted.id == new_rep$predicted.id) / length(orig_rep$rowname)

colnames(orig_rep) <- c("cell_id", "original")
colnames(new_rep) <- c("cell_id", "new")

df <- full_join(orig_rep, new_rep)
df %>% mutate(same = ifelse(original == new, "same", 'different')) %>%
        group_by(original, same) %>%
        summarise(total = n()) %>%
        mutate(percent = total / sum(total)) %>%
        filter(same == "same") %>%
        ggplot(aes(x = original,y =  percent, fill = original)) +
        geom_bar(stat = 'identity') + 
        scale_fill_manual(values = colors)


## Combined UMAP plot

pbmc2$celltype_mrna <- pbmc2_new$predicted.id
pbmc2$celltype_repair <- pbmc2_repair_new$predicted.id
pbmc2$celltype_both <- pbmc2_comb_new$predicted.id

pbmc2 <- RunUMAP(pbmc2, reduction = 'pca', dims = 1:10, spread = 2, min.dist = 1.5)

plot_grid(DimPlot(pbmc2, reduction = 'umap', group.by = 'celltype', cols = colors) + 
                  ggtitle("Reference cell types"),
          DimPlot(pbmc2, reduction = 'umap', group.by = 'celltype_mrna', cols = colors) + 
                  ggtitle("mRNA only reference - 93% the same"),
          DimPlot(pbmc2, reduction = 'umap', group.by = 'celltype_repair', cols = colors) + 
                  ggtitle("Repair only reference - 43% the same"),
          DimPlot(pbmc2, reduction = 'umap', group.by = 'celltype_both', cols = colors) +
                  ggtitle("mRNA and repair reference - 93% the same"))

ggsave("plots/pbmc_classification_plots/all_umaps_simple_celltypes.pdf", 
       height = 8, width = 10, units = 'in', useDingbats = F)

df <- rownames_to_column(FetchData(pbmc2, vars = c("celltype", 
                                                   "celltype_mrna", 
                                                   "celltype_repair", 
                                                   "celltype_both")))

df %>% mutate(mrna = ifelse(celltype == celltype_mrna, 'same', 'diff'),
              repair = ifelse(celltype == celltype_repair, 'same', 'diff'),
              both = ifelse(celltype == celltype_both, 'same', 'diff')) %>%
        select(celltype, mrna, repair, both) %>%
        gather(metric, result, -celltype) %>%
        group_by(celltype, metric, result) %>%
        summarise(total = n()) %>% 
        mutate(percent = total / sum(total)) %>%
        filter(result == "same") -> df

df %>% ggplot(aes(x = celltype, y = percent, fill = celltype)) + 
        geom_bar(stat = 'identity') + 
        theme_cowplot() + 
        scale_fill_manual(values = colors) + 
        facet_grid(~metric) + 
        ylab("Percent of cells identified the same") + 
        xlab("Cell type by reference data") + 
        theme(legend.position = 'none')

ggsave("plots/pbmc_classification_plots/percent_the_same.pdf",
       height = 3, width = 10, units = 'in', useDingbats = F)


### Adding seurat clusters and adding expected data

set.seed(42)
df <- rownames_to_column(FetchData(pbmc2, vars = c("celltype",
                                                   "seurat_clusters",
                                                   "celltype_mrna", 
                                                   "celltype_repair", 
                                                   "celltype_both"))) %>%
        mutate(random_celltype = sample(celltype, ncol(pbmc2)))

cluster_rename = tribble(~seurat_clusters, ~celltype_cluster,
                         0, "NK", 
                         1, "Mono", 
                         2, "T", 
                         3, "T", 
                         4, "T", 
                         5, "B", 
                         6, "Mono", 
                         7, "DC", 
                         8, "Platelet",
                         9, "T", 
                         10, "T") %>%
        mutate(seurat_clusters = as.factor(seurat_clusters))

df %>% left_join(cluster_rename) %>%
        mutate(mrna = ifelse(celltype == celltype_mrna, 'same', 'diff'),
               repair = ifelse(celltype == celltype_repair, 'same', 'diff'),
               both = ifelse(celltype == celltype_both, 'same', 'diff'),
               clusters = ifelse(celltype == celltype_cluster, 'same', 'diff'),
               random = ifelse(celltype == random_celltype, 'same', 'diff'),
               ref = 'same') %>%
        select(celltype, mrna, repair, both, clusters, random, ref) %>%
        gather(metric, result, -celltype) %>%
        group_by(celltype, metric, result) %>%
        summarise(total = n()) %>% 
        mutate(percent = total / sum(total)) %>%
        filter(result == "same") -> df2

c <- c(colors[1:4],colors[6],colors[5])


bar <- function(df) {
        ggplot(df, aes(x = celltype, y = percent, fill = celltype)) + 
                geom_bar(stat = 'identity') + 
                theme_cowplot() + 
                scale_fill_manual(values = c) + 
                ylab("Proportion of cells identified the same") + 
                xlab("Cell type by reference data") + 
                theme(legend.position = 'none') + 
                ylim(0, 1)
}

df2 %>% bar(.) + 
        facet_grid(~metric)

ggsave("plots/pbmc_classification_plots/percent_the_same_withclustersand_random.pdf",
       height = 3, width = 19, units = 'in', useDingbats = F)


df2 %>% group_by(metric) %>%
        nest() %>%
        mutate(plot = map2(data, metric, ~ bar(.x) + ggtitle(.y))) -> plots

plot_grid(plotlist =  plots$plot)

ggsave("plots/pbmc_classification_plots/percent_the_same_withclustersand_random_sep.pdf",
       height = 7, width = 13, units = 'in', useDingbats = F)

# get percentage of sameness
# 93% the same
sum(df$celltype == df$celltype_cluster) / nrow(df)

# 31% the same
sum(df$celltype == df$random_celltype) / nrow(df)

df %>% left_join(cluster_rename) -> df
df <- as.data.frame(df)
rownames(df) <- df$rowname

pbmc2$celltype_cluster <- df$celltype_cluster
pbmc2$celltype_random <- df$random_celltype

plot_grid(DimPlot(pbmc2, reduction = 'umap', group.by = 'celltype', cols = c) + 
                  ggtitle("Reference cell types"),
          DimPlot(pbmc2, reduction = 'umap', group.by = 'celltype_cluster', cols = c) + 
                  ggtitle("Renaming clusters to cell types - 93 % the same"),
          DimPlot(pbmc2, reduction = 'umap', group.by = 'celltype_mrna', cols = c) + 
                  ggtitle("mRNA only reference - 93% the same"),
          DimPlot(pbmc2, reduction = 'umap', group.by = 'celltype_repair', cols = c) + 
                  ggtitle("Repair only reference - 43% the same"),
          DimPlot(pbmc2, reduction = 'umap', group.by = 'celltype_both', cols = c) +
                  ggtitle("mRNA and repair reference - 93% the same"),
          DimPlot(pbmc2, reduction = 'umap', group.by = 'celltype_random', cols = c) + 
                  ggtitle("Random cell types - 31% the same"))

ggsave("plots/pbmc_classification_plots/all_umaps_simple_celltypes.pdf", 
       height = 7, width = 13, units = 'in', useDingbats = F)


df %>% group_by(celltype) %>%
        summarise(number = n()) %>%
        ungroup() %>%
        mutate(percent =  number / sum(number))
df %>% group_by(celltype_repair) %>%
        summarise(number = n()) %>%
        ungroup() %>%
        mutate(percent =  number / sum(number))


## Making confidence intervals for random cell types

ct <- rownames_to_column(FetchData(pbmc2, vars = c("celltype"))) 
n = nrow(ct)

get_percent <- function(df, n, test){
        df <- mutate(df, new = sample(celltype, n),
                     result = ifelse(celltype == new, 'same', 'diff'))
        total <- df %>% group_by(result) %>%
                summarise(n_cells = n()) %>%
                mutate(percent = n_cells/sum(n_cells)) %>%
                filter(result == 'same') %>% 
                mutate(celltype = "all")
        
        df <- df %>% mutate(result = ifelse(celltype == new, 'same', 'diff')) %>%
                select(celltype, result) %>%
                group_by(celltype, result) %>%
                summarise(n_cells = n()) %>%
                mutate(percent = n_cells/sum(n_cells)) %>%
                filter(result == "same") %>% 
                full_join(total) %>% 
                mutate(test = test)
        
        df
}

set.seed(42)
res <- map(c(1:1000), ~get_percent(ct, n, .x))

res <- reduce(res, rbind)
stats <- group_by(res, celltype) %>%
        summarise(average = mean(percent), 
                  stdev = sd(percent),
                  n = 1000) %>%
        mutate(error = qnorm(0.975)*stdev/sqrt(n),
               left = average - error,
               right = average + error
        )
stats %>% filter(celltype != "all") -> stat_df

stat_df <- stat_df %>% mutate(percent = average) 
stat_df %>% bar() + 
        geom_errorbar(aes(ymin=left, ymax=right), width=.2,
                      position=position_dodge(.9))  + 
        ylim(0, 1) + 
        ggtitle("random") -> p1


plot_grid(plots$plot[[1]],
          plots$plot[[2]],
          plots$plot[[3]],
          p1,
          plots$plot[[5]],
          plots$plot[[6]])

ggsave("plots/pbmc_classification_plots/percent_the_same_withclustersand_random_sep.pdf",
       height = 7, width = 13, units = 'in', useDingbats = F)

# Calculating positive and negatives
set.seed(42)
df <- rownames_to_column(FetchData(pbmc2, vars = c("celltype",
                                                   "seurat_clusters",
                                                   "celltype_mrna", 
                                                   "celltype_repair", 
                                                   "celltype_both"))) %>%
        mutate(random_celltype = sample(celltype, ncol(pbmc2)))

cluster_rename = tribble(~seurat_clusters, ~celltype_cluster,
                         0, "NK", 
                         1, "Mono", 
                         2, "T", 
                         3, "T", 
                         4, "T", 
                         5, "B", 
                         6, "Mono", 
                         7, "DC", 
                         8, "Platelet",
                         9, "T", 
                         10, "T") %>%
        mutate(seurat_clusters = as.factor(seurat_clusters))

df %>% left_join(cluster_rename) %>%
        select(-seurat_clusters) %>% 
        gather(method, new_celltype, -rowname, -celltype) -> df

# Get true positive rates e.g. mono -> mono
df %>% mutate(true_positive = ifelse(celltype == new_celltype, T, F)) %>%
        group_by(celltype, method) %>%
        summarise(true_pos = sum(true_positive)) -> true_pos

# Get false positive rates t -> mono
df %>% mutate(neg = ifelse(celltype != new_celltype, T, F)) %>%
        group_by(new_celltype, method) %>%
        summarise(false_pos = sum(neg)) -> false_pos
colnames(false_pos) <- c("celltype", "method", "false_pos")

# Get false negative rates mono -> t
df %>% mutate(neg = ifelse(celltype != new_celltype, T, F)) %>%
        group_by(celltype, method) %>%
        summarise(false_neg = sum(neg)) -> false_neg

# Get true negative rate for mono: anything but mono  > anything but mono
df %>%  mutate(neg = ifelse(celltype != "Mono" & new_celltype != "Mono", T, F)) %>% 
        group_by(method) %>%
        summarise(true_neg = sum(neg)) %>%
        mutate(celltype = "Mono") -> true_neg_mono

# Get true negative rate for t: anything but t  > anything but t
df %>%  mutate(neg = ifelse(celltype != "T" & new_celltype != "T", T, F)) %>% 
        group_by(method) %>%
        summarise(true_neg = sum(neg)) %>%
        mutate(celltype = "T") -> true_neg_T

# Get true negative rate for nk: anything but nk  > anything but nk
df %>%  mutate(neg = ifelse(celltype != "NK" & new_celltype != "NK", T, F)) %>% 
        group_by(method) %>%
        summarise(true_neg = sum(neg)) %>%
        mutate(celltype = "NK") -> true_neg_NK

# Get true negative rate for dc: anything but dc  > anything but dc
df %>%  mutate(neg = ifelse(celltype != "DC" & new_celltype != "DC", T, F)) %>% 
        group_by(method) %>%
        summarise(true_neg = sum(neg)) %>%
        mutate(celltype = "DC") -> true_neg_DC

# Get true negative rate for platelet: anything but platelet  > anything but platelet
df %>% mutate(neg = ifelse(celltype != "Platelet" & new_celltype != "Platelet", T, F)) %>% 
        group_by(method) %>%
        summarise(true_neg = sum(neg)) %>%
        mutate(celltype = "Platelet") -> true_neg_P

# Get true negative rate for B: anything but B  > anything but B
df %>%  mutate(neg = ifelse(celltype != "B" & new_celltype != "B", T, F)) %>% 
        group_by(method) %>%
        summarise(true_neg = sum(neg)) %>%
        mutate(celltype = "B") -> true_neg_B

# Combining true negative rates
true_neg <- reduce(list(true_neg_DC, true_neg_mono, true_neg_NK,
                        true_neg_P, true_neg_T, true_neg_B), rbind)

# Combining all results
test_res <- reduce(list(true_pos, false_neg, false_pos, true_neg), full_join)

# Calculating precision, recall, true positive and true neg rates
test_res <- test_res %>%
        mutate(precision = true_pos / (true_pos + false_pos),
               recall = true_pos / (true_pos + false_neg),
               TP_rate = true_pos / (true_pos + false_neg),
               FP_rate = false_pos / (false_pos + true_neg))

# Plot precision vs recall
test_res %>% 
        ggplot(aes(x = precision, y = recall, color = celltype)) +
        geom_point(aes(shape = method))  + 
        scale_color_manual(values = c)

# Plot true positive vs false positive
test_res %>% 
        ggplot(aes(x = FP_rate, y = TP_rate, color = celltype)) +
        geom_point(aes(shape = method),size = 5)  + 
        scale_color_manual(values = c) + 
        xlim(0 ,1) + 
        geom_abline(slope = 1, intercept = 0, color = "#999999") + 
        theme_minimal_grid() + 
        ylab("True Positive Rate") + 
        xlab("False Positive Rate")

ggsave("plots/pbmc_classification_plots/ROCish_plot.pdf",
       height = 4, width = 6, units = 'in', useDingbats = F)

# Plot positive predictive value for eacc method
test_res %>%
        ggplot(aes(x = method, y = precision, fill = celltype)) + 
        geom_bar(stat = 'identity', position = 'dodge') + 
        scale_fill_manual(values = c) + 
        theme_cowplot()

ggsave("plots/pbmc_classification_plots/PPV_barplot.pdf",
       height = 3, width = 10, units = 'in', useDingbats = F)




