## Supplementary Figure 3

source("scripts/functions.R")
library(tidyverse) 
library(scrunchy)
library(Seurat)

load("../data/barnyard/barnyard.seurat.object.Rdata")

DimPlot(object = barnyard_seurat, reduction = 'umap', group.by = "seurat_clusters",
        cols = colors)

ggsave(filename = "../plots/barnyard/barnyard_UMAP_seurat.pdf", height = 4, width = 4, units = 'in',
       useDingbats = F)


## Read in colors determined in Fig1

color_df <- read_tsv("../data/barnyard/celltypes_by_repair.tsv")
color_df <- as.data.frame(color_df)
rownames(color_df) <- color_df$cell_id

barnyard_seurat$celltype <- color_df["color"]

df = rownames_to_column(FetchData(barnyard_seurat, 
                                  vars = c("celltype", "seurat_clusters", 
                                                            "UNG", "RNASEH2C",
                                                            "Uracil2-45", "riboG2-44")),
                        "cell_id")
colnames(df) <- str_replace(colnames(df), "-", "_")

# Plot expression and repair activity by cluster
# plot_activity is a function of scrunchy

p1 <- plot_activity(df, UNG, seurat_clusters) + 
        theme(legend.position = "none") + 
        ggtitle("UNG mRNA expression") + 
        labs(x = "mRNA expression")
p2 <- plot_activity(df, repair_Uracil2_45, seurat_clusters) + 
        theme(legend.position = "none") + 
        ggtitle("Uracil repair activity")
p3 <- plot_activity(df, RNASEH2C, seurat_clusters) + 
        theme(legend.position = "none") + 
        ggtitle("RNASEH2C mRNA expression") + 
        labs(x = "mRNA expression")
p4 <- plot_activity(df, repair_riboG2_44, seurat_clusters) + 
        theme(legend.position = "none") + 
        ggtitle("Ribonucleotide repair activity")

plot_grid(p1, p3, p2, p4, nrow = 1)

ggsave("../plots/barnyard/expression_and_activity_beeswarm_by_cluster.pdf",
       height = 6, width = 7.5, units = 'in',useDingbats = F)


# Make table of cells expression mRNA or have repair activity

df %>% mutate(UNG_expression = ifelse(UNG == 0, 
                                      "0", "> 0"),
              RNASEH2C_expression = ifelse(RNASEH2C == 0,
                                           "0", "> 0"),
              Uracil_repair = ifelse(repair_Uracil2_45 == 0, 
                                     "0", "> 0"),
              Ribonucleotide_repair = ifelse(repair_riboG2_44 == 0,
                                             "0", "> 0")) %>%
        select(seurat_clusters, UNG_expression, RNASEH2C_expression, 
               Uracil_repair, Ribonucleotide_repair) %>%
        gather(measurement, status, -seurat_clusters) %>%
        group_by(seurat_clusters, measurement, status) %>%
        summarise(count = n()) %>%
        ungroup() %>%
        filter(status != "0") %>%
        select(-status) %>%
        mutate(measurement = str_replace(measurement, "_", " ")) %>%
        spread(measurement, count) %>%
        select(seurat_clusters, `UNG expression`, `RNASEH2C expression`, 
               `Uracil repair`, `Ribonucleotide repair`) -> number_cells

total_cells <- df %>% group_by(seurat_clusters) %>%
        summarise(`Total number of cells` = n())

table <- full_join(number_cells, total_cells)
colnames(table) <- c("Cluster",
                     "UNG expression",
                     "RNASEH2C expression",
                     "Uracil repair",
                     "Ribonucleotide repair",
                     "Total number of cells")

# Save tables as pdf using gt package
library(gt)
gt1 <- gt(table)
gtsave(gt1, filename = "../plots/barnyard/expression_and_repair_cellnumbers.pdf")

table %>% gather(measurment, count, -Cluster) %>%
        group_by(measurment) %>%
        summarise(Total = sum(count)) %>%
        spread(measurment, Total) %>%
        mutate(Cluster = 'Total') %>% gt() -> gt2

gtsave(gt2, filename = "../plots/barnyard/expression_and_repair_total_numbers.pdf")

## Making barnyard plots colored by cluster and cell types

# Rename clusters based on repair of that cluster
cluster_rename = tribble(~seurat_clusters, ~celltype_from_cluster,
                         0, "RNASEH2CKO", 
                         1, "RNASEH2CKO", 
                         2, "UNGKO", 
                         3, "UNGKO",
                         4, "UNGKO", 
                         5, "RNASEH2CKO",
                         6, "Both") %>%
        mutate(seurat_clusters = as.factor(seurat_clusters))

df %>% full_join(cluster_rename) -> df

# Getting count data from seurat object for plotting barnyard
count_df <- rownames_to_column(FetchData(barnyard_seurat, 
                                         vars = c("Uracil2-45", "riboG2-44", "seurat_clusters"),
                                         slot = 'counts'),
                               "cell_id")

df %>% select(cell_id, celltype, celltype_from_cluster) %>%
        full_join(count_df) -> count_df


count_df %>% ggplot(aes(x = `repair_Uracil2-45`, 
                        y = `repair_riboG2-44`, 
                        color = celltype_from_cluster)) + 
        ggbeeswarm::geom_quasirandom(groupOnX = F,alpha = 0.7) + 
        scale_color_manual(values = colors[1:3]) + 
        theme_cowplot()  + 
        theme(legend.position = "top",
              legend.title = element_blank()) + 
        ggtitle("Celltype determined from clusters") + 
        xlab("Counts at uracil repair site") + 
        ylab("Counts at ribonucleotide repair site") -> p2

count_df %>% ggplot(aes(x = `repair_Uracil2-45`, 
                        y = `repair_riboG2-44`, 
                        color = seurat_clusters)) + 
        ggbeeswarm::geom_quasirandom(groupOnX = F, alpha = 0.7) + 
        scale_color_manual(values = colors) + 
        theme_cowplot()  + 
        theme(legend.position = "top",
              legend.title = element_blank()) + 
        xlab("Counts at uracil repair site") + 
        ylab("Counts at ribonucleotide repair site") + 
        ggtitle("Clusters") -> p1

count_df %>% ggplot(aes(x = `repair_Uracil2-45`, 
                        y = `repair_riboG2-44`, 
                        color = celltype)) + 
        ggbeeswarm::geom_quasirandom(groupOnX = F, alpha = 0.7) + 
        scale_color_manual(values = c(colors[1], "#999999", colors[2:3])) + 
        theme_cowplot()  + 
        theme(legend.position = "top",
              legend.title = element_blank()) + 
        xlab("Counts at uracil repair site") + 
        ylab("Counts at ribonucleotide repair site") + 
        ggtitle("Celltype determined by repair counts")-> p3

plot_grid(p1, p2, p3, nrow = 1)

ggsave("../plots/barnyard/barnyard_plot_dodge_clusters_andcelltypes.pdf",
       height = 4, width = 12, units = 'in', useDingbats = F)




