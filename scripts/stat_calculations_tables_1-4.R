## Calculating statistics for repair and gene expression differences

library(Seurat)
library(tidyverse)
library(clustifyr)
library(cowplot)
source("scripts/functions.R")

load("../data/pbmc/seurat/pbmc1.seurat.Rdata")
pbmc1 <- subset(pbmc1, subset = celltype != "Platelet")
load("../data/pbmc/seurat/pbmc2.seurat.Rdata")
pbmc2 <- subset(pbmc2, subset = celltype != "Platelet")
load("../data/pbmc/seurat/pbmc3.seurat.Rdata")
pbmc3 <- subset(pbmc3, subset = celltype != "Platelet")


## Calculating cell type vs all other cells for RNA and repair
Idents(pbmc1) <- pbmc1$celltype
mark1 <- rownames_to_column(FindAllMarkers(pbmc1, assay = "RNA"), "cell_id")
dna1 <- rownames_to_column(FindAllMarkers(pbmc1, assay = "repair"), "cell_id")

Idents(pbmc2) <- pbmc2$celltype
mark2 <- rownames_to_column(FindAllMarkers(pbmc2, assay = "RNA"), "cell_id")
dna2 <- rownames_to_column(FindAllMarkers(pbmc2, assay = "repair"), "cell_id")

Idents(pbmc3) <- pbmc3$celltype
mark3 <- rownames_to_column(FindAllMarkers(pbmc3, assay = "RNA"), "cell_id")
dna3 <- rownames_to_column(FindAllMarkers(pbmc3, assay = "repair"), "cell_id")

## Combing results into one table
marks <- purrr::map2(list(mark1, mark2, mark3), c("pbmc1", "pbmc2", "pbmc3"),
                     ~ mutate(.x, sample = .y))
dna <- purrr::map2(list(dna1, dna2, dna3), c("pbmc1", "pbmc2", "pbmc3"),
                   ~ mutate(.x, sample = .y))

marks <- reduce(marks, rbind)
dna <- reduce(dna, rbind)


## Filter for base exicion repair genes from KEGG database: KEGG id ko03410
repair_genes <- read_tsv("../data/other/BER_genes.txt") %>% 
        filter(KEGG_BASE_EXCISION_REPAIR != "> Base excision repair" )
repair_genes <- repair_genes$KEGG_BASE_EXCISION_REPAIR
repair_genes <- c(repair_genes, "RNASEH2A", "RNASEH2B", "RNASEH2C")

# filtering for significanct BER genes only
marks %>% filter(gene %in% repair_genes) -> ber

# filtering for repair positions only
repair_pos <- c("Uracil-45", "GU-45", "Abasic-45", "Abasic-46", "riboG-44")
dna %>% filter(gene %in% repair_pos,
               p_val_adj < 0.05) -> rep_dna

# plots for stats
plot_bar <- function(df, t) {
        ggplot(df, aes(x = gene, y = avg_logFC)) + 
                geom_bar(stat = "identity") + 
                geom_text(aes(label = signif(p_val_adj, 2))) + 
                theme_cowplot() + 
                ylab("Average log fold change") + 
                geom_hline(yintercept = 0, alpha= .5) + 
                xlab(NULL) + 
                theme(legend.position = 'none') + 
                scale_fill_manual(values = colors) + 
                ggtitle(t)
}


ber %>% ggplot(aes(x = gene, y = avg_logFC)) + 
        geom_bar(stat = 'identity') + 
        facet_grid(sample~cluster)

## Dendritic cells only
marks %>% filter(sample == "pbmc1", cluster == "DC",
                 gene %in% repair_genes) %>%
        plot_bar(. , "pbmc1") -> p1
marks %>% filter(sample == "pbmc2", cluster == "DC",
                 gene %in% repair_genes) %>%
        plot_bar(., "pbmc2") -> p2
marks %>% filter(sample == "pbmc3", cluster == "DC",
                 gene %in% repair_genes) %>%
        plot_bar(., "pbmc3") -> p3
plot_grid(p1, p2, p3, nrow = 1)

# Split list to make all plots
ber %>% group_by(cluster, sample) %>%
        group_split() -> berList

# Plot all sample and celltype stats separately 
map2(berList, unique(paste(ber$cluster, ber$sample)), ~plot_bar(.x, .y)) -> pl 
plot_grid(plotlist = pl)

## Merge seurat object and find significant differences between cell types for all cells
all <- merge(pbmc1, c(pbmc2, pbmc3))
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 5000)  
all <- ScaleData(all)
all <- RunPCA(all, npcs = 20, verbose = FALSE)
all <- RunUMAP(all, reduction = "pca", dims = 1:20)
all <- FindNeighbors(all, reduction = "pca", dims = 1:20)
all <- FindClusters(all, resolution = 0.5)

DimPlot(all, reduction = 'umap', group.by = "celltype")
Idents(all) <- all$celltype

# Find differentially expressed RNA
all_mark <- FindAllMarkers(all, assay = "RNA")

# Find differential repair 
all_dna <- FindAllMarkers(all, assay = "repair")

# Filter for repair positions
all_dna %>% filter(gene %in% repair_pos, 
                   p_val_adj < 0.05) -> all_dna_rep

# Filter for repair genes
all_mark %>% filter(gene %in% repair_genes) -> all_rep

#Split for plotting
all_rep %>% group_by(cluster) %>%
        group_split -> allList

# Make all by all plots
map2(allList, unique(all_rep$cluster), 
     ~plot_bar(.x, .y)) -> pl 

plot_grid(plotlist = pl)

# Save tsv files for RNA expression repair tables are saved below
write_tsv(all_mark, path = "../data/tables/all_pbmc_replicates_markers.tsv")
write_tsv(marks, path = "../data/tables/all_pbmc_replicates_separate_markers.tsv")
write_tsv(all_rep, path = "../data/tables/all_pbmc_replicates_markers_repair.tsv")
write_tsv(ber, path = "../data/tables/all_pbmc_replicates_separate_markers_repair.tsv")

## Calculating pairwise differences for REPAIR stats: 

# Make list of pairwise differences
comb = t(combn(unique(pbmc1$celltype), 2))
comb = data.frame(comb[, 1], comb[,2])
colnames(comb) <- c("var1", "var2")
comb <- mutate(comb, var1 = as.character(var1),
               var2 = as.character(var2))

# Funcrtion to calculate stat and return df
mark <- function(seurat, var1, var2, ...){
        res = rownames_to_column(FindMarkers(seurat, 
                                             ident.1 = var1, ident.2 = var2, ...),
                                 'gene')
        res = mutate(res, celltype1 = var1, celltype2 = var2)
        res
}

# Calculate stats for all pairwise comparisons
res = map2(comb$var1, comb$var2, ~ mark(pbmc1, .x, .y, assay = 'repair'))
res = reduce(res, rbind)

# Filter of significant adjusted p values
res %>% filter(p_val_adj < 0.05) %>%
        mutate(sample = "pbmc1") -> df

# Repeat for other samples
res = map2(comb$var1, comb$var2, ~ mark(pbmc2, .x, .y, assay = 'repair'))
res = reduce(res, rbind)

res %>% filter(p_val_adj < 0.05) %>%
        mutate(sample = "pbmc2") %>%
        full_join(df) -> df

res = map2(comb$var1, comb$var2, ~ mark(pbmc3, .x, .y, assay = 'repair'))
res = reduce(res, rbind)

res %>% filter(p_val_adj < 0.05) %>%
        mutate(sample = "pbmc3") %>%
        full_join(df) -> df

# Repeat for all cells together
res = map2(comb$var1, comb$var2, ~ mark(all, .x, .y, assay = 'repair'))
res = reduce(res, rbind)

res %>% filter(p_val_adj < 0.05) %>%
        mutate(sample = "all_cells") %>%
        full_join(df) -> df

# Filter for repair positions and save
df %>% filter(gene %in% repair_pos) %>%  
        arrange(sample, gene, p_val_adj) %>%
        write_tsv(., path = "../data/tables/all_pbmcs_REPAIR_data_only.tsv")

# Filter for hairpins presenet in all samples (AU, GU, riboG, Abasic, Normal)
# TI and CI were not included
df %>% separate(gene, into = c("hairpin", "position")) %>%
        filter(position != 'biotin',
               hairpin %in% c("Uracil", "GU", 
                              "riboG", "Abasic", "Normal")) %>%
        unite(gene, hairpin, position, sep = "-") %>%
        arrange(sample, gene, p_val_adj) %>%
        write_tsv(., path = "../data/tables/all_pbmcs_REPAIR_all_pos.tsv")

## Filtering celltype vs all cells for hairpins present in all samples and 
## Combing tables made above 
all_dna %>% mutate(sample = "all_cells") %>%
        full_join(dna) %>% 
        select(-cell_id) %>%
        separate(gene, into = c("hairpin", "position")) %>%
        filter(position != 'biotin',
               hairpin %in% c("Uracil", "GU", 
                              "riboG", "Abasic", "Normal"),
               p_val_adj < 0.05) %>%
        unite(gene, hairpin, position, sep = "-") %>%
        select( "gene",
                "p_val",
                "avg_logFC",
                "pct.1",
                "pct.2",
                "p_val_adj",
                "cluster",
                "sample" ) %>%
        arrange(sample, gene, p_val_adj) %>%
        write_tsv(., path = "../data/tables/all_pbmcs_REPAIR_all_pos_vsallcells.tsv")

## Filtering celltype vs all cells for hairpins present in all samples
## Combing tables made above 

all_dna %>% mutate(sample = "all_cells") %>%
        full_join(dna) %>% 
        select(-cell_id) %>%
        filter(gene %in% repair_pos) %>%
        select( "gene",
                "p_val",
                "avg_logFC",
                "pct.1",
                "pct.2",
                "p_val_adj",
                "cluster",
                "sample" ) %>%
        arrange(sample, gene, p_val_adj) %>%
       write_tsv(., path = "../data/tables/all_pbmcs_REPAIR_repair_pos_vsallcells.tsv")



