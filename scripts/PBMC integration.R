## Integrating of PBMC samples and creating figures

source("scripts/functions.R")

load("../data/pbmc/seurat/pbmc1.seurat.Rdata")
pbmc1$sample <- 1
load("../data/pbmc/seurat/pbmc2.seurat.Rdata")
pbmc2$sample <- 2
load("../data/pbmc/seurat/pbmc3.seurat.Rdata")
pbmc3$sample <- 3

immune.anchors <- FindIntegrationAnchors(object.list = c(pbmc1, pbmc2, pbmc3), dims = 1:20)
pbmc.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

#Run PCA
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)
pbmc.combined <- RunPCA(pbmc.combined, npcs = 20, verbose = FALSE)
# t-SNE and Clustering
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca", dims = 1:20)
pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca", dims = 1:20)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.5)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = 'pca', dims = 1:10, spread = 2, min.dist = 1.5)

## Do the same for the repair integrated data
repair.anchors <- FindIntegrationAnchors(object.list = c(pbmc1, pbmc2, pbmc3), 
                                         assay = c("repair", "repair", "repair"), 
                                         dims = 1:20)

pbmc.repair <- IntegrateData(anchorset = repair.anchors, dims = 1:20)

#Run PCA
pbmc.repair <- ScaleData(pbmc.repair, verbose = FALSE)
pbmc.repair <- RunPCA(pbmc.repair, npcs = 20, verbose = FALSE)
# t-SNE and Clustering
pbmc.repair <- RunUMAP(pbmc.repair, reduction = "pca", dims = 1:20)
pbmc.repair <- FindNeighbors(pbmc.repair, reduction = "pca", dims = 1:20)
pbmc.repair <- FindClusters(pbmc.repair, resolution = 0.5)
pbmc.repair <- RunUMAP(pbmc.repair, reduction = 'pca', dims = 1:10, spread = 2, min.dist = 1.5)

# Should integrate based on mRNA
plot_grid(DimPlot(pbmc.combined, reduction = 'umap', group.by = "celltype"),
          DimPlot(pbmc.repair, reduction = 'umap', group.by = "celltype"),
          DimPlot(pbmc.combined, reduction = 'umap', group.by = "sample"),
          DimPlot(pbmc.repair, reduction = 'umap', group.by = "sample"))

pbmc.combined = NormalizeData(object = pbmc.combined, assay = "repair", method = "LogNormalize")
pbmc.combined = ScaleData(object = pbmc.combined, assay = "repair")
save(pbmc.combined, file = "../data/pbmc/seurat/pbmc_123_integrated.seurat.Rdata")


pbmc.combined <- subset(pbmc.combined, subset = celltype != "Platelet")

plot_grid(DimPlot(pbmc.combined, reduction = 'umap', group.by = "celltype", cols = colors),
          DimPlot(pbmc.combined, reduction = 'umap', group.by = "sample", cols = colors))
          
ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/integrated_UMAP.pdf",
       height = 4, width = 8, units = 'in', useDingbats = F)

## Making bulk coverage plots
df <- get_hairpin_coverage(pbmc.combined)

df %>% mutate(adduct_position1 = 44,
              adduct_position2 = -1) %>%
        filter(celltype != "Platelet") -> df

df %>%
        filter(hairpin == 'Uracil',
               position > 34) %>%
        haircut_plot(., x = "position", y = "avg_count", point = TRUE,
                     xlim = c(35, 55), pal = colors, col = 'celltype', 
                     y_lab = "Average counts per cell") + 
        theme(legend.position = 'top') + 
        ggtitle("U:A repair")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/integrated_UA_bulk.pdf",
       height = 4, width = 5, units = 'in', useDingbats = F)

df %>%
        filter(hairpin == 'GU',
               position > 34) %>%
        haircut_plot(., x = "position", y = "avg_count", point = TRUE,
                     xlim = c(35, 55), pal = colors, col = 'celltype', 
                     y_lab = "Average counts per cell") + 
        theme(legend.position = 'top') + 
        ggtitle("U:G repair")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/integrated_GU_bulk.pdf",
       height = 4, width = 5, units = 'in', useDingbats = F)


df %>%
        filter(hairpin == 'riboG',
               position > 34) %>%
        haircut_plot(., x = "position", y = "avg_count", point = TRUE,
                     xlim = c(35, 55), pal = colors, col = 'celltype', 
                     y_lab = "Average counts per cell") + 
        theme(legend.position = 'top') + 
        ggtitle("Ribonucleotide repair")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/integrated_riboG_bulk.pdf",
       height = 4, width = 5, units = 'in', useDingbats = F)


df %>%
        filter(hairpin == 'Abasic',
               position > 34) %>%
        haircut_plot(., x = "position", y = "avg_count", point = TRUE,
                     xlim = c(35, 55), pal = colors, col = 'celltype', 
                     y_lab = "Average counts per cell") + 
        theme(legend.position = 'top') + 
        ggtitle("Abasic repair")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/integrated_Ab_bulk.pdf",
       height = 4, width = 5, units = 'in', useDingbats = F)


df %>%
        filter(hairpin == 'Normal',
               position > 34) %>%
        haircut_plot(., x = "position", y = "avg_count", point = TRUE,
                     xlim = c(35, 55), pal = colors, col = 'celltype', 
                     y_lab = "Average counts per cell") + 
        theme(legend.position = 'top') + 
        ggtitle("Unmodified substrate")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/integrated_Normal_bulk.pdf",
       height = 4, width = 5, units = 'in', useDingbats = F)

## Single cell plots
repair.positions = c("Uracil-45", 
                     "riboG-44", 
                     "GU-45", 
                     "Abasic-46", 
                     "Abasic-45", 
                     "Normal-45")
df <- get_single_cell_df(pbmc.combined, feat = c(repair.positions, "celltype"))

#Make tidy data
df %>% gather(repair, activity, -celltype, -cell_id) -> df

#Add labels for plotting
repair_labels = tribble(~repair, ~label,
                        "Uracil_45", "U:A repair",
                        "GU_45", "U:G repair",
                        "riboG_44", "Ribonucelotide repair",
                        "Abasic_46", "Abasic repair long-patch",
                        "Abasic_45", "Abasic repair short-patch",
                        "Normal_45", "Unmodified substrate"
)

# Put samples in correct order
df %>% full_join(repair_labels) %>%        
        mutate(label = fct_relevel(label, "U:A repair", 
                                   "U:G repair", 
                                   "Ribonucelotide repair", 
                                   "Abasic repair long-patch",
                                   "Abasic repair short-patch",
                                   "Unmodified substrate")) -> df

df %>% filter(repair %in% c("Uracil_45", "GU_45", "riboG_44", "Normal_45")) %>%
        activity_plot() + 
        facet_wrap(~label, ncol = 1, strip.position = "left") 

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/integrated_single_cell.pdf",
       height = 8, width = 4, units = "in", useDingbats = F)

df %>% filter(repair %in% c("Abasic_45", "Abasic_46")) %>%
        activity_plot(lab = label) 

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/integrated_single_cell_ab.pdf",
       height = 3, width = 4, units = "in", useDingbats = F)

## Calculating stats

Idents(pbmc.combined) <- pbmc.combined$celltype

rna <- rownames_to_column(FindAllMarkers(pbmc.combined, assay = "RNA"), "cell_id")
repair <- rownames_to_column(FindAllMarkers(pbmc.combined, assay = "repair"), "cell_id")

# Filter for repair genes
repair_genes <- read_tsv("../data/other/BER_genes.txt") %>% 
        filter(KEGG_BASE_EXCISION_REPAIR != "> Base excision repair" )
repair_genes <- repair_genes$KEGG_BASE_EXCISION_REPAIR
repair_genes <- c(repair_genes, "RNASEH2A", "RNASEH2B", "RNASEH2C")

# filtering for significanct BER genes only
rna %>% filter(gene %in% repair_genes) -> ber

# filtering for significant repair sites
repair_pos <- c("Uracil-45", "GU-45", "Abasic-45", "Abasic-46", "riboG-44")
repair %>% filter(gene %in% repair_pos,
               p_val_adj < 0.05) -> repair_only

write_tsv(rna, path = "../data/tables/integrated_pbmc_mrna_markers.tsv")
write_tsv(ber, path = "../data/tables/integrated_pbmc_mrna_ber_markers.tsv")
write_tsv(repair, path = "../data/tables/integrated_pbmc_repair_markers.tsv")
write_tsv(repair_only, path = "../data/tables/integrated_pbmc_repair_sitesonly_markers.tsv")


## Calculating pairwise differences for REPAIR stats: 
comb = t(combn(unique(pbmc.combined$celltype), 2))
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
res = map2(comb$var1, comb$var2, ~ mark(pbmc.combined, .x, .y, assay = 'repair'))
res = reduce(res, rbind)

# Filter of significant adjusted p values
res %>% filter(p_val_adj < 0.05) -> df

df %>% filter(gene %in% repair_pos) %>%  
        mutate(sample = "integrated") %>%
        arrange(gene, p_val_adj) %>%
        write_tsv("../data/tables/integrated_pbmc_pairwise_repair.tsv")
# Filter for hairpins presenet in all samples (AU, GU, riboG, Abasic, Normal)
# TI and CI were not included
df %>% separate(gene, into = c("hairpin", "position")) %>%
        mutate(sample = "integrated") %>%
        filter(position != 'biotin',
               hairpin %in% c("Uracil", "GU", 
                              "riboG", "Abasic", "Normal")) %>%
        unite(gene, hairpin, position, sep = "-") %>%
        arrange(sample, gene, p_val_adj) %>%
        write_tsv(., path = "../data/tables/integrated_pbmc_pairwise_all_pos.tsv")








