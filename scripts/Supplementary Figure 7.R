## Supplementary Figure 7
## Calculate number of oligos captured per drop compared to number of oligos in drop

source("scripts/functions.R")

## Barnyard data

load("../data/barnyard/barnyard.seurat.object.Rdata")

repair_position = tribble(~hairpin, ~repair_position, ~amt_added, 
                          ~adduct_position1, ~adduct_position2, ~number_added,
                          "Uracil1", 45, 10, 44, 0, 800000,
                          "Uracil2", 45, 25, 44, 0, 2000000,
                          "Uracil3", 45, 5, 44, 0, 400000,
                          "Uracil4", 45, 2.5, 44, 0, 200000,
                          "Uracil5", 45, 0.5, 44, 0, 40000,
                          "riboG1", 44, 10, 44, 0, 800000,
                          "riboG2", 44, 25, 44, 0, 2000000,
                          "riboG3", 44, 5,44, 0, 400000,
                          "riboG4", 44, 2.5,44, 0, 200000,
                          "riboG5", 44, 0.5, 44, 0, 40000,
                          "CI", 45, 5, 44, 0, 400000,
                          "TI", 45, 5, 44, 0, 400000,
                          "Abasic", 45, 5, 44, 0, 400000,
                          "Normal", 1, 5, 44, 0, 400000,
                          "GU", 45, 5, 44, 0, 400000)

get_repair_df <- function(object, assay = 'repair', ...){
        m <- t(as.matrix(GetAssayData(object, assay = assay, ...)))
        df <- rownames_to_column(as.data.frame(m), "cell_id")
        df <- gather(df, hp_pos, value, -cell_id) %>%
                separate(hp_pos, into = c("hairpin", "position"), sep = '-')
        df
}

color_df <- read_tsv("../data/barnyard/celltypes_by_repair.tsv")
color_df <- as.data.frame(color_df)
rownames(color_df) <- color_df$cell_id

barnyard_seurat$celltype <- color_df["color"]

dat <- get_repair_df(barnyard_seurat, slot = 'counts')
dat <- dat %>% left_join(repair_position)

dat %>% group_by(hairpin, cell_id) %>%
        mutate(total_count = sum(value)) -> dat

dat %>% filter(position == repair_position) %>%
        select(cell_id, hairpin, position, number_added, total_count, value, amt_added) %>% 
        left_join(color_df) -> dat

dat %>% mutate(captured = total_count/number_added,
               repaired = value / total_count,
               repaired = replace(repaired, is.na(repaired), 0))-> prop

prop %>% ungroup() %>% mutate(hairpin = str_remove(hairpin, '[0-9]+')) %>%
        filter(amt_added == 5,
               color %in% c("RNASEH2KO", "UNGKO")) %>%
        ggplot(aes(x = hairpin, y = captured, color = hairpin)) + 
        ggbeeswarm::geom_quasirandom() + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        ylab("Proportion of hairpins captured") + 
        theme(legend.position = 'none') + 
        xlab(NULL) -> p1

ggsave(plot = p1, filename = "~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/barnyard_prop_of_hairpins_captured_all_hp_allcells.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)


prop %>% ungroup() %>% mutate(hairpin = str_remove(hairpin, '[0-9]+')) %>%
        filter(amt_added == 5,
               color %in% c("RNASEH2KO", "UNGKO")) %>%
        ggplot(aes(x = hairpin, y = captured, color = hairpin)) + 
        ggbeeswarm::geom_quasirandom() + 
        facet_grid(~color) + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        ylab("Proportion of hairpins captured") + 
        theme(legend.position = 'none') -> p2
ggsave(plot = p2, filename = "~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/barnyard_prop_of_hairpins_captured_all_hp_allcells_split.pdf",
       height = 4, width = 8, units = 'in', useDingbats = F)


prop %>% ungroup() %>% mutate(hairpin = str_remove(hairpin, '[0-9]+')) %>%
        filter(amt_added == 5,
               color %in% c("RNASEH2KO", "UNGKO")) %>%
        ggplot(aes(x = hairpin, y = repaired, color = hairpin)) + 
        ggbeeswarm::geom_quasirandom() + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        ylab("Proportion of captured hairipns repaired") + 
        theme(legend.position = 'none') -> p3

ggsave(plot = p3, filename = "~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/barnyard_prop_of_hairpins_repaired_all_hp_allcells.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)


prop %>% ungroup() %>% mutate(hairpin = str_remove(hairpin, '[0-9]+')) %>%
        filter(amt_added == 5,
               color %in% c("RNASEH2KO", "UNGKO")) %>%
        ggplot(aes(x = hairpin, y = repaired, color = hairpin)) + 
        ggbeeswarm::geom_quasirandom() + 
        facet_grid(~color) + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        ylab("Proportion of captured hairipns repaired") + 
        theme(legend.position = 'none') -> p4

ggsave(plot = p4, filename = "~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/barnyard_prop_of_hairpins_repaired_all_hp_allcells_split.pdf",
       height = 4, width = 8, units = 'in', useDingbats = F)

plot_grid(p1, p2, p3, p4, ncol = 2, rel_widths = c(1, 1.5))

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/barnyard_prop_of_hairpins_captured_all_hp.pdf",
       height = 8, width = 16, units = 'in', useDingbats = F)

prop %>% ungroup() %>% mutate(hairpin = str_remove(hairpin, '[0-9]+')) %>%
        filter(hairpin %in% c("Uracil", "riboG"),
               color %in% c("RNASEH2KO", "UNGKO")) %>%
        ggplot(aes(x = as.factor(amt_added), y = captured, color = hairpin)) + 
        ggbeeswarm::geom_quasirandom(dodge.width = 1) + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        ylab("Proportion of captured hairpins") + 
        xlab("Concentration of hairpin in master mix (nM)") + 
        theme(legend.position = "top",
              legend.title = element_blank()) -> p1

ggsave(plot = p1, filename = "~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/barnyard_prop_of_hairpins_captured_dilution_allcells.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)


prop %>% ungroup() %>% mutate(hairpin = str_remove(hairpin, '[0-9]+')) %>%
        filter(hairpin %in% c("Uracil", "riboG"),
               color %in% c("RNASEH2KO", "UNGKO")) %>%
        ggplot(aes(x = as.factor(amt_added), y = captured, color = hairpin)) + 
        ggbeeswarm::geom_quasirandom(dodge.width = 1) + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        facet_grid(~color) + 
        ylab("Proportion of captured hairpins") + 
        xlab("Concentration of hairpin in master mix (nM)") + 
        theme(legend.position = "top",
              legend.title = element_blank()) -> p2

ggsave(plot = p2, filename = "~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/barnyard_prop_of_hairpins_captured_dilution_allcells_splt.pdf",
       height = 4, width = 8, units = 'in', useDingbats = F)

prop %>% ungroup() %>% mutate(hairpin = str_remove(hairpin, '[0-9]+')) %>%
        filter(hairpin %in% c("Uracil", "riboG"),
               color %in% c("RNASEH2KO", "UNGKO")) %>%
        ggplot(aes(x = as.factor(amt_added), y = repaired, color = hairpin)) + 
        ggbeeswarm::geom_quasirandom(dodge.width = 1) + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        ylab("Proportion of captured hairpins repaired") + 
        xlab("Concentration of hairpin in master mix (nM)") + 
        theme(legend.position = "top",
              legend.title = element_blank()) -> p3

ggsave(plot = p3, filename = "~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/barnyard_prop_of_hairpins_repaired_dilution_allcells.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)

prop %>% ungroup() %>% mutate(hairpin = str_remove(hairpin, '[0-9]+')) %>%
        filter(hairpin %in% c("Uracil", "riboG"),
               color %in% c("RNASEH2KO", "UNGKO")) %>%
        ggplot(aes(x = as.factor(amt_added), y = repaired, color = hairpin)) + 
        ggbeeswarm::geom_quasirandom(dodge.width = 1) + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        facet_grid(~color) + 
        ylab("Proportion of captured hairpins repaired") + 
        xlab("Concentration of hairpin in master mix (nM)") + 
        theme(legend.position = "top",
              legend.title = element_blank()) -> p4

ggsave(plot = p4, filename = "~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/barnyard_prop_of_hairpins_repaired_dilution_allcells_split.pdf",
       height = 4, width = 8, units = 'in', useDingbats = F)

plot_grid(p1, p2, p3, p4, ncol = 2, rel_widths = c(1, 1.5))

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/barnyard_prop_of_hairpins_captured_dilution.pdf",
       height = 8, width = 16, units = 'in', useDingbats = F)

## Calculating statstics

prop %>% filter(amt_added == 5,
                color %in% c("RNASEH2KO", "UNGKO")) %>%
        select(hairpin, color, captured, cell_id) %>%
        spread(hairpin, captured) %>%
        ungroup() %>%
        select(-cell_id) -> stat_5

stat_activity_grouped(stat_5, color) %>% View()

prop %>% select(cell_id, hairpin, captured, repaired) %>%
        gather(prop, value, -cell_id, -hairpin) %>%
        unite(measurement, hairpin, prop) %>%
        spread(measurement, value) -> meta

meta <- as.data.frame(meta) 
rownames(meta) <- meta$cell_id
meta <- select(meta, -cell_id)
meta<- as.matrix(meta)
meta <- t(meta)

barnyard_seurat[["proportion"]] <- CreateAssayObject(counts = meta)
Idents(barnyard_seurat) <- barnyard_seurat$celltype

prop %>% ungroup() %>% select(captured, repaired, hairpin) -> stat_df
stat_activity_grouped(stat_df, hairpin) %>% View()

r <- tribble(~hairpin, ~status,
             "Uracil3", "repaired", 
             "Abasic", "repaired", 
             "riboG3", "repaired", 
             "GU", "repaired", 
             "CI", "not_repaired",
             "TI", "not_repaired",
             "Normal", "not_repaired")

prop %>% ungroup() %>% filter(amt_added == 5) %>%
        select(captured, repaired, hairpin) %>% 
        left_join(r) %>%
        select(-hairpin) -> stat_df
stat_activity_grouped(stat_df, status) %>% View()


# No significant differences between groups (celltypes)
FindAllMarkers(barnyard_seurat, assay = "proportion")




## PBMCs

load("../data/pbmc/seurat/pbmc60.seurat.Rdata")
load("../data/pbmc/seurat/pbmc30.seurat.Rdata")
load("../data/pbmc/seurat/pbmc15.seurat.Rdata")

s <- tribble(~sample,
               "pbmc15", 
               "pbmc30",
               "pbmc60")
        
dat <- mutate(s,
              contents = map(c(pbmc15, pbmc30, pbmc60), 
                             ~get_repair_df(.x, slot = 'counts')))

dat <- unnest(dat)
dat <- dat %>% left_join(repair_position)

ct <- mutate(s, 
             cellt = map(c(pbmc15, pbmc30, pbmc60),
                           ~ get_single_cell_df(.x, feat = c("celltype")) %>%
                                 mutate(cell_id = str_remove(cell_id, "_1")))) %>%
        unnest()


dat <- dat %>% left_join(ct)

dat %>% group_by(sample, hairpin, cell_id) %>%
        mutate(total_count = sum(value)) -> dat

dat %>% filter(position == repair_position) %>%
        select(sample, cell_id, hairpin, position, 
               number_added, total_count, 
               value, amt_added, celltype) -> dat

reads <- tribble(~sample, ~read_num,
                 'pbmc15', 13568589,
                 'pbmc30', 13334806,
                 'pbmc60', 13289873)

dat %>% select(sample, cell_id) %>%
        group_by(sample) %>%
        distinct(cell_id) %>%
        summarise(number_of_cells = n()) -> num_of_cells

dat %>% left_join(reads) %>% 
        left_join(num_of_cells) %>%
        mutate(reads_per_cell = read_num / number_of_cells) -> dat

dat %>% group_by(sample, cell_id) %>%
        mutate(total_algined_reads = sum(value),
               prop_aligned = total_algined_reads / reads_per_cell) %>%
        select(sample, cell_id, total_algined_reads, prop_aligned) %>%
        distinct() -> aligned_reads
        
aligned_reads %>% ggplot(aes(x = sample, y = prop_aligned)) + 
        ggbeeswarm::geom_quasirandom() + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        ylab("Proportion of aligned reads per cell") + 
        xlab("sample") + 
        theme(legend.position = "top",
              legend.title = element_blank())

dat %>% mutate(captured = total_count/number_added,
               repaired = value / total_count,
               repaired = replace(repaired, is.na(repaired), 0))-> prop

prop %>% ungroup() %>% mutate(hairpin = str_remove(hairpin, '[0-9]+')) %>%
        filter(hairpin %in% c("Uracil", "riboG")) %>% 
        separate(sample, into = c("sample", "Time"), sep = -2) -> dil

dil %>% ggplot(aes(x = as.factor(amt_added), y = captured, color = Time)) + 
        ggbeeswarm::geom_quasirandom(dodge.width = 1) + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        ylab("Proportion of captured hairpins") + 
        xlab("Concentration of hairpin in master mix (nM)") + 
        theme(legend.position = "top") + 
        facet_grid(~hairpin) 

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/pbmc_dilution_time_proportion_captured.pdf",
       height = 4, width = 8, units = 'in', useDingbats = F)

dil %>% ggplot(aes(x = as.factor(amt_added), y = repaired, color = Time)) + 
        ggbeeswarm::geom_quasirandom(dodge.width = 1) + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        ylab("Proportion of captured hairpins repaired") + 
        xlab("Concentration of hairpin in master mix (nM)") + 
        theme(legend.position = "top") + 
        facet_grid(~hairpin)

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/pbmc_dilution_time_proportion_repaired.pdf",
       height = 4, width = 8, units = 'in', useDingbats = F)


dil %>% ggplot(aes(x = as.factor(amt_added), y = captured, color = as.factor(time))) + 
        ggbeeswarm::geom_quasirandom(dodge.width = 1) + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        ylab("Proportion of captured hairpins") + 
        xlab("Concentration of hairpin in master mix (nM)") + 
        theme(legend.position = "top",
              legend.title = element_blank()) + 
        facet_grid(hairpin~celltype)

dil %>% ggplot(aes(x = as.factor(amt_added), y = repaired, color = as.factor(time))) + 
        ggbeeswarm::geom_quasirandom(dodge.width = 1) + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        ylab("Proportion of captured hairpins repaired") + 
        xlab("Concentration of hairpin in master mix (nM)") + 
        theme(legend.position = "top",
              legend.title = element_blank()) + 
        facet_grid(hairpin~celltype)


ggsave(plot = p1, filename = "~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/pbmc_prop_of_hairpins_captured_dilution_hp_allcells.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)


prop %>% ungroup() %>% mutate(hairpin = str_remove(hairpin, '[0-9]+')) %>%
        filter(hairpin %in% c("Uracil", "riboG")) %>%
        ggplot(aes(x = as.factor(amt_added), y = captured, color = hairpin)) + 
        ggbeeswarm::geom_quasirandom(dodge.width = 1) + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        facet_grid(~celltype) + 
        ylab("Proportion of captured hairpins") + 
        xlab("Concentration of hairpin in master mix (nM)") + 
        theme(legend.position = "top",
              legend.title = element_blank()) -> p2

ggsave(plot = p2, filename = "~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/pbmc_prop_of_hairpins_captured_dilution_allcells_split.pdf",
       height = 4, width = 8, units = 'in', useDingbats = F)


prop %>% ungroup() %>% mutate(hairpin = str_remove(hairpin, '[0-9]+')) %>%
        filter(hairpin %in% c("Uracil", "riboG")) %>%
        ggplot(aes(x = as.factor(amt_added), y = repaired, color = hairpin)) + 
        ggbeeswarm::geom_quasirandom(dodge.width = 1) + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        ylab("Proportion of captured hairpins repaired") + 
        xlab("Concentration of hairpin in master mix (nM)") + 
        theme(legend.position = "top",
              legend.title = element_blank()) -> p3

ggsave(plot = p3, filename = "~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/pbmc_prop_of_hairpins_repaired_dilution_hp_allcells.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)


prop %>% ungroup() %>% mutate(hairpin = str_remove(hairpin, '[0-9]+')) %>%
        filter(hairpin %in% c("Uracil", "riboG")) %>%
        ggplot(aes(x = as.factor(amt_added), y = repaired, color = hairpin)) + 
        ggbeeswarm::geom_quasirandom(dodge.width = 1) + 
        theme_cowplot() + 
        scale_color_manual(values = colors) + 
        facet_grid(~celltype) + 
        ylab("Proportion of captured hairpins repaired") + 
        xlab("Concentration of hairpin in master mix (nM)") + 
        theme(legend.position = "top",
              legend.title = element_blank()) -> p4

ggsave(plot = p4, filename = "~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/pbmc_prop_of_hairpins_repaired_dilution_split.pdf",
       height = 4, width = 8, units = 'in', useDingbats = F)


plot_grid(p1, p2, p3, p4, ncol = 2, rel_widths = c(1, 4))

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/pbmc_prop_of_hairpins_captured_dilution.pdf",
       height = 8, width = 24, units = 'in', useDingbats = F)


