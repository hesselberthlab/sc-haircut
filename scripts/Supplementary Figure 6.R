# Supplementary Figure 6
## Dilution and timecourse of PBMC data
source("scripts/functions.R")

# Load seurat objects
load("../data/pbmc/seurat/pbmc15.seurat.Rdata")
load("../data/pbmc/seurat/pbmc30.seurat.Rdata")
load("../data/pbmc/seurat/pbmc60.seurat.Rdata")

# Add repair positions and amount added data
repair_position = tribble(~hairpin, ~repair_position, ~amt_added, ~adduct_position1, ~adduct_position2,
                          "Uracil1", 45, 10, 44, 0,
                          "Uracil2", 45, 25, 44, 0,
                          "Uracil3", 45, 5, 44, 0,
                          "Uracil4", 45, 2.5, 44, 0,
                          "Uracil5", 45, 0.5, 44, 0,
                          "riboG1", 44, 10, 44, 0,
                          "riboG2", 44, 25, 44, 0,
                          "riboG3", 44, 5,44, 0, 
                          "riboG4", 44, 2.5,44, 0,
                          "riboG5", 44, 0.5, 44, 0,
                          "CI", 45, 5, 44, 0,
                          "TI", 45, 5, 44, 0,
                          "Abasic", 45, 5, 44, 0,
                          "Normal", 1, 5, 44, 0,
                          "GU", 45, 5, 44, 0)

get_repair_df <- function(object, assay = 'repair', ...){
        m <- t(as.matrix(GetAssayData(object, assay = assay, ...)))
        df <- rownames_to_column(as.data.frame(m), "cell_id")
        df <- gather(df, hp_pos, value, -cell_id) %>%
                separate(hp_pos, into = c("hairpin", "position"), sep = '-')
        df
}

# Get celltypes
celltypes <- map(c(pbmc15, pbmc30, pbmc60),
                 ~get_single_cell_df(.x, feat = c("celltype", "seurat_clusters")))
celltypes <- map(celltypes, ~mutate(.x, cell_id = str_remove(cell_id, "_1")))

# Make dataframe
norm_df <- map(c(pbmc15, pbmc30, pbmc60), get_repair_df)
norm_df <- map2(norm_df, c(15, 30, 60), ~mutate(.x, time = .y))
norm_df <- map2(norm_df, celltypes, left_join)
norm_df <- reduce(norm_df, full_join)
norm_df <- norm_df %>% left_join(repair_position)

colnames(norm_df) <- str_remove(colnames(norm_df), "_seurat")

#Filter for repair positions
repair_df <- norm_df %>% filter(position == repair_position)

# Get diltuion dataframe
dldf <- repair_df %>% 
        filter(!hairpin %in% c("CI", "TI", "Normal", "GU", "Abasic")) %>%
        separate(hairpin, into = c("hairpin", "number"), sep = -1) 

#Make plots
dldf %>% filter(hairpin == "Uracil",
                celltype != "Platelet") %>%
        ggplot(aes(y = as.factor(amt_added), x = value, color = celltype)) + 
        ggbeeswarm::geom_quasirandom(groupOnX = F) + 
        facet_grid(celltype~time) + 
        scale_color_manual(values = colors) + 
        cowplot::theme_cowplot() + 
        theme(legend.position = 'none',
              strip.background.y = element_blank()) + 
        geom_vline(xintercept = c(2, 4, 6), alpha = 0.5, 
                   color = "#999999", linetype = "dashed") + 
        xlab("Activity") + 
        ylab("Substrate concentration (nM)")

#Save plot
ggsave("~/hesselberthlab/projects/10x_haircut/20190819/plots/pbmc/uracil_repair_beswarm_concentration_time_newcelltype.pdf",
       height = 8, width = 9, units = 'in', useDingbats = F)

#Make plot
dldf %>% filter(hairpin == "riboG",
                celltype != "Platelet") %>%
        ggplot(aes(y = as.factor(amt_added), x = value, color = celltype)) + 
        ggbeeswarm::geom_quasirandom(groupOnX = F) + 
        facet_grid(celltype~time) + 
        scale_color_manual(values = colors) + 
        cowplot::theme_cowplot() + 
        theme(legend.position = 'none',
              strip.background.y = element_blank()) + 
        geom_vline(xintercept = c(2, 4, 6), alpha = 0.5, 
                   color = "#999999", linetype = "dashed") + 
        xlab("Activity") + 
        ylab("Substrate concentration (nM)")

#Save plot
ggsave("~/hesselberthlab/projects/10x_haircut/20190819/plots/pbmc/riboG_repair_beswarm_concentration_time_newcelltype.pdf",
       height = 8, width = 9, units = 'in', useDingbats = F)

dldf %>% filter(celltype != "Platelet",
                amt_added == 25) %>%
        ggplot(aes(y = as.factor(time), x = value, color = celltype)) + 
        ggbeeswarm::geom_quasirandom(groupOnX = F) + 
        facet_grid(celltype~hairpin) + 
        scale_color_manual(values = colors) + 
        cowplot::theme_cowplot() + 
        theme(legend.position = 'none',
              strip.background.y = element_blank()) + 
        geom_vline(xintercept = c(2, 4, 6), alpha = 0.5, 
                   color = "#999999", linetype = "dashed") + 
        xlab("Activity") + 
        ylab("Substrate concentration (nM)")


ggsave("~/hesselberthlab/projects/10x_haircut/20190819/plots/pbmc/uracil_riboG_25nm_bytime.pdf",
       height = 8, width = 6, units = 'in', useDingbats = F)

#Save plot
ggsave("~/hesselberthlab/projects/10x_haircut/20190819/plots/pbmc/uracil_repair_beswarm_concentration_time_newcelltype.pdf",
       height = 8, width = 9, units = 'in', useDingbats = F)


### Slope plot

dldf %>% select(-number, -position, 
                -adduct_position1, 
                -adduct_position2,
                -repair_position) %>%
        spread(hairpin, value) -> dldf_spread

# Funciton to calculate slope
lm_f <- function(df){
        res <- lm( Uracil ~ amt_added, df)
        res <- tidy(res)
        res <- filter(res, term == 'amt_added')
        res <- mutate(res, substrate = "uracil")
        id <- df[["cell_id"]][1]
        res <- mutate(res, cell_id = id)
        res
}

# Calculate uracil slope
dldf_spread %>% group_by(time, cell_id, celltype) %>%
        group_modify(~ broom::tidy(lm(Uracil ~ amt_added, data = .x))) %>%
        mutate(substrate = "Uracil") -> u_stat

#Calculate ribo slope
dldf_spread %>% group_by(time, cell_id, celltype) %>%
        group_modify(~ broom::tidy(lm(riboG ~ amt_added, data = .x))) %>%
        mutate(substrate = "riboG")-> r_stat

# Make df
stat_df <- full_join(u_stat, r_stat)

#Make plot
stat_df %>% filter(term == "amt_added",
                   celltype != "Platelet") %>% 
        ggplot(aes(x = estimate, y = as.factor(time), color = celltype)) + 
        ggbeeswarm::geom_quasirandom(groupOnX = FALSE) + 
        scale_color_manual(values = colors) + 
        facet_grid(celltype~substrate) + 
        cowplot::theme_cowplot() + 
        theme(legend.position = 'none',
              strip.background.y = element_blank()) + 
        ylab("Time (min)") + 
        xlab("Slope of normalized amount vs amt added") + 
        geom_vline(xintercept = c(0.0, 0.1, 0.2), alpha = 0.5, linetype = 'dashed',color = "#999999")


#Save plot
#ggsave("~/hesselberthlab/projects/10x_haircut/20190819/plots/pbmc/slope_beeswarm_bytime_newcelltypes.pdf",
#       height = 8, width = 6, units = 'in', useDingbats = F)

## Add linear model to seurat object and test significance 
stat_df %>% filter(term == "amt_added", 
                   time == 15) %>% 
        ungroup() %>%
        select(cell_id, estimate, substrate) %>%
        spread(substrate, estimate) -> slope15
        
slope15 <- as.data.frame(slope15)
rownames(slope15) <- slope15$cell_id
slope15 <- slope15[, 2:3 ]
slope15 <- as.matrix(slope15)

pbmc15$slope <- CreateAssayObject(t(slope15))

stat_df %>% filter(term == "amt_added", 
                   time == 30) %>% 
        ungroup() %>%
        select(cell_id, estimate, substrate) %>%
        spread(substrate, estimate) -> slope30

slope30 <- as.data.frame(slope30)
rownames(slope30) <- slope30$cell_id
slope30 <- slope30[, 2:3 ]
slope30 <- as.matrix(slope30)

pbmc30$slope <- CreateAssayObject(t(slope30))

stat_df %>% filter(term == "amt_added", 
                   time == 60) %>% 
        ungroup() %>%
        select(cell_id, estimate, substrate) %>%
        spread(substrate, estimate) -> slope60

slope60 <- as.data.frame(slope60)
rownames(slope60) <- slope60$cell_id
slope60 <- slope60[, 2:3 ]
slope60 <- as.matrix(slope60)

pbmc60$slope <- CreateAssayObject(t(slope60))

Idents(pbmc15) <- pbmc15$celltype
pbmc15 <- subset(pbmc15, subset = celltype != "Platelet")
Idents(pbmc30) <- pbmc30$celltype
pbmc30 <- subset(pbmc30, subset = celltype != "Platelet")
Idents(pbmc60) <- pbmc60$celltype
pbmc60 <- subset(pbmc60, subset = celltype != "Platelet")

dna <- purrr::map2(list(pbmc15, pbmc30, pbmc60), c(15, 30, 60), ~rownames_to_column(
        FindAllMarkers(.x, assay = "slope", logfc.threshold = 0), "cell_id") %>%
                mutate(time = .y))

dna <- reduce(dna, rbind)
## Not significant genes this way

comb = t(combn(unique(pbmc15$celltype), 2))
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

## Calculate pairwise for each time
res = map2(comb$var1, comb$var2, ~ mark(pbmc15, .x, .y, assay = 'slope',
                                        logfc.threshold = 0))
res = reduce(res, rbind)

# Filter of significant adjusted p values
res %>% filter(p_val_adj < 0.05) %>%
        mutate(time = 15) -> df

# Repeat for other samples
res = map2(comb$var1, comb$var2, ~ mark(pbmc30, .x, .y, assay = 'slope',
                                        logfc.threshold = 0))
res = reduce(res, rbind)

res %>% filter(p_val_adj < 0.05) %>%
        mutate(time = 30) %>%
        full_join(df) -> df

res = map2(comb$var1, comb$var2, ~ mark(pbmc60, .x, .y, assay = 'slope',
                                        logfc.threshold = 0))
res = reduce(res, rbind)

res %>% filter(p_val_adj < 0.05) %>%
        mutate(time = 60) %>%
        full_join(df) -> df

df %>% arrange(gene, time, p_val_adj) %>% View()
        write_tsv("../data/tables/timecourse_pbmc_slope_pairwise.tsv")

dna %>% arrange(gene, time, p_val_adj) %>% 
        write_tsv("../data/tables/timecourse_pbmc_Slope_allothercells.tsv")


# Make and save umap plots
pbmc15 <- subset(pbmc15, subset = celltype != "Platelet")
DimPlot(pbmc15, reduction = 'umap', cols = colors, group.by = "celltype")
#ggsave("plots/pbmc/UMAP_15min_new_celltype.pdf",
#       height = 4, width = 5, units = 'in', useDingbats = F)

pbmc30 <- subset(pbmc30, subset = celltype != "Platelet")
DimPlot(pbmc30, reduction = 'umap', cols = colors, group.by = "celltype")
#ggsave("~/hesselberthlab/projects/10x_haircut/20190819/plots/pbmc/UMAP_30min_new_celltype.pdf",
#       height = 4, width = 5, units = 'in', useDingbats = F)

pbmc60 <- subset(pbmc60, subset = celltype != "Platelet")
DimPlot(pbmc60, reduction = 'umap', cols = colors, group.by = "celltype")
#ggsave("plots/pbmc/UMAP_60min_new_celltype.pdf",
#       height = 4, width = 5, units = 'in', useDingbats = F)
