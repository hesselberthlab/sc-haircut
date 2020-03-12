# Supplementary Figure 3
## Bulk plots vs empty drops -- finding biological signal

library(tidyverse)
library(Seurat)
library(cowplot)
library(scrunchy)
source("scripts/functions.R")

# Read in cell data
cells <- Read10X("../data/pbmc/pbmc1/filtered_pbmc_haircut/",
                gene.column = 1)

# Read in empty drops from repair pipeline (before filtering)
empty <- Read10X("../data/pbmc/pbmc1/pbmc1_haircut/", gene.column = 1)

# Filter out drops with cells
empty <- empty[, !colnames(empty) %in% colnames(cells)]

# Make matrices into bulk dataframes
celldf <- as.data.frame(rownames(cells))
celldf$counts <- Matrix::rowSums(cells)
celldf$drop <- rep("cell", nrow(celldf))
colnames(celldf) <- c("hp_pos", "count", "drop")

# Make df total counts per drop for empty drops
df <- as.data.frame(colnames(empty))
df$counts <- Matrix::colSums(empty)

# Sort df by total counts per cell and get the top ~4000 drops
df %>% arrange(desc(counts)) %>% slice(1:ncol(cells)) -> empty4000

# Make empty drop matrix into bulk dataframe
emptydf <- as.data.frame(rownames(empty))
emptydf$counts <- Matrix::rowSums(empty[, empty4000$`colnames(empty)`])
emptydf$drop <- rep("empty", nrow(emptydf))
colnames(emptydf) <- c("hp_pos", "count", "drop")

# Combind df
df <- rbind(celldf, emptydf)

n <- ncol(cells)

# Get average counts per position
df %>% separate(hp_pos, into = c("hairpin", "position")) %>%
        mutate(position = as.double(position),
               avg_count = count /n) -> df

# Make bulk plots
c = c(colors[1], "#999999")
haircut_plot(df, x = "position", y = "avg_count", col = "drop", 
             point = T, pal = c) + 
        facet_grid(~hairpin)

# Make bulk plots - separately 
df %>% group_by(hairpin) %>% nest() %>%
        mutate(plot = map2(data, hairpin, ~haircut_plot(.x, x = "position", 
                                       y = "avg_count", 
                                       col = "drop", 
                                       point = T, pal = c, 
                                       y_lab = "Average counts per drop") + 
                                   ggtitle(.y))) -> plots

plot_grid(plotlist = plots$plot)

#Save plots
#ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/bulk_emptydrops_allsub/bulk_emptydrops_pbmc1.pdf",
#       height = 9, width = 10, units = 'in', useDingbats = F)

## Calculating stats for repair vs not repair

# Get matrix of empy ~4000 drops
e <- empty[, empty4000$`colnames(empty)`]
number <- ncol(e)

#Turn empty and cell matrices into dataframes
e <- rownames_to_column(as.data.frame(e), "hp_pos") %>%
        gather(cellid, count, -hp_pos) %>%
        mutate(drop = 'empty')
sc_df <- rownames_to_column(as.data.frame(cells), "hp_pos") %>%
        gather(cellid, count, -hp_pos) %>%
        mutate(drop = 'cell') %>%
        full_join(e)

#Function to get wilcoson results
tidy_wilcoxon <- function(x, y) {
        broom::tidy(suppressWarnings(stats::wilcox.test(x, y)))
}

# get stats per position comparing cells to empty drops by position
sc_df %>% arrange(drop) %>%
        mutate(number = rep(1:1707834, 2)) %>%
        select(-cellid) %>%
        spread(drop, count) %>%
        group_by(hp_pos) %>%
        group_modify(~ tidy_wilcoxon(.x$cell, .x$empty)) -> pval

# Calculate qvalues
 qval <- qvalue::qvalue(pval$p.value, pi0 = 1)
 
 # add qvalues to pval results
 pval %>% ungroup() %>%
         mutate(q.value = qval$qvalues) -> pval
 
 pval %>% filter(qvalue < 0.05)
 # All sites are significant qvalues except for 31 sites

## PBMC2 bulk plots
cells <- Read10X("../data/pbmc/pbmc2/filtered_pbmc_haircut/",
                 gene.column = 1)

empty <- Read10X("../data/pbmc/pbmc2/pbmc2_haircit/", gene.column = 1)
empty <- empty[, !colnames(empty) %in% colnames(cells)]

celldf <- as.data.frame(rownames(cells))
celldf$counts <- Matrix::rowSums(cells)
celldf$drop <- rep("cell", nrow(celldf))
colnames(celldf) <- c("hp_pos", "count", "drop")


df <- as.data.frame(colnames(empty))
df$counts <- Matrix::colSums(empty)

df %>% arrange(desc(counts)) %>% slice(1:ncol(cells)) -> empty4000

emptydf <- as.data.frame(rownames(empty))
emptydf$counts <- Matrix::rowSums(empty[, empty4000$`colnames(empty)`])
emptydf$drop <- rep("empty", nrow(emptydf))
colnames(emptydf) <- c("hp_pos", "count", "drop")

df <- rbind(celldf, emptydf)

n <- ncol(cells)

df %>% mutate(hp_pos = str_replace(hp_pos, "_biotin", "Biotin")) %>%
        separate(hp_pos, into = c("hairpin", "position")) %>%
        mutate(position = as.double(position),
               avg_count = count /n) -> df

c = c(colors[1], "#999999")
haircut_plot(df, x = "position", y = "avg_count", col = "drop", 
             point = T, pal = c) + 
        facet_grid(~hairpin)

df %>% group_by(hairpin) %>% nest() %>%
        mutate(plot = map2(data, hairpin, ~haircut_plot(.x, x = "position", 
                                                        y = "avg_count", 
                                                        col = "drop", 
                                                        point = T, pal = c, 
                                                        y_lab = "Average counts per drop") + 
                                   ggtitle(.y))) -> plots

plot_grid(plotlist = plots$plot)
#ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/bulk_emptydrops_allsub/bulk_emptydrops_pbmc2.pdf",
#       height = 9, width = 13, units = 'in', useDingbats = F)

# PBMC3
cells <- Read10X("../data/pbmc/pbmc3/filtered_pbmc_haircut/",
                 gene.column = 1)

empty <- Read10X("../data/pbmc/pbmc3/pbmc3_haircut/", gene.column = 1)
empty <- empty[, !colnames(empty) %in% colnames(cells)]

celldf <- as.data.frame(rownames(cells))
celldf$counts <- Matrix::rowSums(cells)
celldf$drop <- rep("cell", nrow(celldf))
colnames(celldf) <- c("hp_pos", "count", "drop")


df <- as.data.frame(colnames(empty))
df$counts <- Matrix::colSums(empty)

df %>% arrange(desc(counts)) %>% slice(1:ncol(cells)) -> empty4000

emptydf <- as.data.frame(rownames(empty))
emptydf$counts <- Matrix::rowSums(empty[, empty4000$`colnames(empty)`])
emptydf$drop <- rep("empty", nrow(emptydf))
colnames(emptydf) <- c("hp_pos", "count", "drop")

df <- rbind(celldf, emptydf)

n <- ncol(cells)

df %>% mutate(hp_pos = str_replace(hp_pos, "_biotin", "Biotin")) %>%
        separate(hp_pos, into = c("hairpin", "position")) %>%
        mutate(position = as.double(position),
               avg_count = count /n) -> df

c = c(colors[1], "#999999")
haircut_plot(df, x = "position", y = "avg_count", col = "drop", 
             point = T, pal = c) + 
        facet_grid(~hairpin)

df %>% group_by(hairpin) %>% nest() %>%
        mutate(plot = map2(data, hairpin, ~haircut_plot(.x, x = "position", 
                                                        y = "avg_count", 
                                                        col = "drop", 
                                                        point = T, pal = c, 
                                                        y_lab = "Average counts per drop") + 
                                   ggtitle(.y))) -> plots

plot_grid(plotlist = plots$plot)
#ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/bulk_emptydrops_allsub/bulk_emptydrops_pbmc3.pdf",
#       height = 9, width = 13, units = 'in', useDingbats = F)


### PST1 sample for MGMT activity

cells <- Read10X("~/hesselberthlab/projects/10x_haircut/20190221/data/haircut/PBM1_PST1/filtered/",
                 gene.column = 1)

empty <- Read10X("~/hesselberthlab/projects/10x_haircut/20190221/data/haircut/PBM1_PST1/", 
                 gene.column = 1)
empty <- empty[, !colnames(empty) %in% colnames(cells)]

celldf <- as.data.frame(rownames(cells))
celldf$counts <- Matrix::rowSums(cells)
celldf$drop <- rep("cell", nrow(celldf))
colnames(celldf) <- c("hp_pos", "count", "drop")


df <- as.data.frame(colnames(empty))
df$counts <- Matrix::colSums(empty)

df %>% arrange(desc(counts)) %>% slice(1:ncol(cells)) -> empty4000

emptydf <- as.data.frame(rownames(empty))
emptydf$counts <- Matrix::rowSums(empty[, empty4000$`colnames(empty)`])
emptydf$drop <- rep("empty", nrow(emptydf))
colnames(emptydf) <- c("hp_pos", "count", "drop")

df <- rbind(celldf, emptydf)

n <- ncol(cells)

df %>% mutate(hp_pos = str_replace(hp_pos, "_biotin", "Biotin"),
              hp_pos = str_replace(hp_pos, "pst1.", '')) %>%
        separate(hp_pos, into = c("hairpin", "position"), sep = "_") %>%
        mutate(position = as.double(position),
               avg_count = count /n) -> df

c = c(colors[1], "#999999")
haircut_plot(df, x = "position", y = "avg_count", col = "drop", 
             point = T, pal = c) + 
        facet_grid(~hairpin)

df %>% group_by(hairpin) %>% nest() %>%
        mutate(plot = map2(data, hairpin, ~haircut_plot(.x, x = "position", 
                                                        y = "avg_count", 
                                                        col = "drop", 
                                                        point = T, pal = c, 
                                                        y_lab = "Average counts per drop") + 
                                   ggtitle(.y))) -> plots

plot_grid(plotlist = plots$plot)
#ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/bulk_emptydrops_allsub/bulk_emptydrops_pbmc2_pst1.pdf",
#       height = 9, width = 13, units = 'in', useDingbats = F)


### PST1 sample for MGMT activity - PBMC3

cells <- Read10X("../data/pbmc/pbmc2_pst1/filtered_haircut",
                 gene.column = 1)

empty <- Read10X("../data/pbmc/pbmc2_pst1/umitools_counts.tsv.gz", 
                 gene.column = 1)
empty <- empty[, !colnames(empty) %in% colnames(cells)]

celldf <- as.data.frame(rownames(cells))
celldf$counts <- Matrix::rowSums(cells)
celldf$drop <- rep("cell", nrow(celldf))
colnames(celldf) <- c("hp_pos", "count", "drop")


df <- as.data.frame(colnames(empty))
df$counts <- Matrix::colSums(empty)

df %>% arrange(desc(counts)) %>% slice(1:ncol(cells)) -> empty4000

emptydf <- as.data.frame(rownames(empty))
emptydf$counts <- Matrix::rowSums(empty[, empty4000$`colnames(empty)`])
emptydf$drop <- rep("empty", nrow(emptydf))
colnames(emptydf) <- c("hp_pos", "count", "drop")

df <- rbind(celldf, emptydf)

n <- ncol(cells)

df %>% mutate(hp_pos = str_replace(hp_pos, "_biotin", "Biotin"),
              hp_pos = str_replace(hp_pos, "pst1.", '')) %>%
        separate(hp_pos, into = c("hairpin", "position"), sep = "_") %>%
        mutate(position = as.double(position),
               avg_count = count /n) -> df

c = c(colors[1], "#999999")
haircut_plot(df, x = "position", y = "avg_count", col = "drop", 
             point = T, pal = c) + 
        facet_grid(~hairpin)

df %>% group_by(hairpin) %>% nest() %>%
        mutate(plot = map2(data, hairpin, ~haircut_plot(.x, x = "position", 
                                                        y = "avg_count", 
                                                        col = "drop", 
                                                        point = T, pal = c, 
                                                        y_lab = "Average counts per drop") + 
                                   ggtitle(.y))) -> plots

plot_grid(plotlist = plots$plot)

#ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/bulk_emptydrops_allsub/bulk_emptydrops_pbmc3_pst1.pdf",
#       height = 9, width = 13, units = 'in', useDingbats = F)


## Reads per drop don't really change with cells
# Read in cell data
cells <- Read10X("../data/pbmc/pbmc1/filtered_pbmc_haircut/",
                 gene.column = 1)

# Read in empty drops from repair pipeline (before filtering)
empty <- Read10X("../data/pbmc/pbmc1/pbmc1_haircut/", gene.column = 1)

per_all <- Matrix::colSums(empty)
# Filter out drops with cells
empty <- empty[, !colnames(empty) %in% colnames(cells)]

per_cell <- Matrix::colSums(cells)
per_drop <- Matrix::colSums(empty)


per_all <- as.data.frame(per_all)
per_all <- arrange(per_all, desc(per_all))
csx <- cumsum(per_all$per_all) 
csx_small = csx[seq(1, length(csx), 10000)]


plot(x = 1:length(csx_small),         # Plot of cumsum vector
     y = csx_small)

#~ 70000 drops with hairpins
csx_small = csx[1:70000]
plot(x = 1:length(csx_small),         # Plot of cumsum vector
     y = csx_small)

true_drops <- per_all[1:70000,]
mean(true_drops)

n_reads <- 54699753
n_cells <- 4009

#reads per cell 13644.24
n_reads / n_cells

# reads per drop 781.425
n_reads / 70000

# Average aligned reads / drop 582.6254
mean(true_drops)

# ~80% of reads align to hairpins 625.14
(n_reads / 70000) * .8

## BARNYARD Data

celltypes <- read_tsv("../data/barnyard/celltypes_by_repair.tsv")
ungko <- filter(celltypes, color == "UNGKO")
ungko <- ungko$cell_id

rnasehko <- filter(celltypes, color == "RNASEH2KO") 
rnasehko <- rnasehko$cell_id

# Read in cell data
cells <- Read10X("../data/barnyard/filtered_haircut/",
                 gene.column = 1)

ungko <- cells[, ungko]
rnasehko <- cells[, rnasehko]

# Read in empty drops from repair pipeline (before filtering)
empty <- Read10X("../data/barnyard/haircut_mtx/", gene.column = 1)

# Filter out drops with cells
empty <- empty[, !colnames(empty) %in% colnames(cells)]

# Make matrices into bulk dataframes
ungkodf <- as.data.frame(rownames(ungko))
ungkodf$counts <- Matrix::rowSums(ungko)
ungkodf$drop <- rep("ungko", nrow(ungko))
colnames(ungkodf) <- c("hp_pos", "count", "drop")


rnasedf <- as.data.frame(rownames(rnasehko))
rnasedf$counts <- Matrix::rowSums(rnasehko)
rnasedf$drop <- rep("rnasehko", nrow(rnasehko))
colnames(rnasedf) <- c("hp_pos", "count", "drop")

celldf <- rbind(ungkodf, rnasedf)

# Make df total counts per drop for empty drops
df <- as.data.frame(colnames(empty))
df$counts <- Matrix::colSums(empty)

# Sort df by total counts per cell and get the top ~4000 drops
df %>% arrange(desc(counts)) %>% slice(1:ncol(cells)) -> empty4000

# Make empty drop matrix into bulk dataframe
emptydf <- as.data.frame(rownames(empty))
emptydf$counts <- Matrix::rowSums(empty[, empty4000$`colnames(empty)`])
emptydf$drop <- rep("empty", nrow(emptydf))
colnames(emptydf) <- c("hp_pos", "count", "drop")

# Combind df
df <- rbind(celldf, emptydf)

n <- tribble(~drop, ~total,
             'ungko', 911,
             'rnasehko', 1307,
             'empty', 2377)

# Get average counts per position
df %>% separate(hp_pos, into = c("hairpin", "position")) %>%
        left_join(n) %>%
        mutate(position = as.double(position),
               avg_count = count /total) -> df

# Make bulk plots
c = c("#999999", colors[c(2, 3)])
haircut_plot(df, x = "position", y = "avg_count", col = "drop", 
             point = T, pal = c) + 
        facet_grid(~hairpin)

# Make bulk plots - separately 
df %>% group_by(hairpin) %>% nest() %>%
        mutate(plot = map2(data, hairpin, ~haircut_plot(.x, x = "position", 
                                                        y = "avg_count", 
                                                        col = "drop", 
                                                        point = T, pal = c, 
                                                        y_lab = "Average counts per drop") + 
                                   ggtitle(.y))) -> plots

plot_grid(plotlist = plots$plot)

ggsave("~/Dropbox (Hesselberth Lab)/HesselberthLab/Images/10x_barnyard_emptry_drops_with_celltypes.pdf",
       height = 10, width = 9, useDingbats = F )



                