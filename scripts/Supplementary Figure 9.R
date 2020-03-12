## Supplementary Figure 9
## DNA repair factor mRNA and repair for ours and reference

source("scripts/functions.R")
load("../data/pbmc/seurat/pbmc1.seurat.Rdata")

## Got reference data from protein atlas
## RNA consensus from http://www.proteinatlas.org
# RNA consensus tissue gene data
# Consensus transcript expression levels summarized per gene 
# in 74 tissues based on transcriptomics data from three sources: 
#         HPA, GTEx and FANTOM5. 
# The consensus normalized expression ("NX") value is calculated as the 
# maximum NX value for each gene in the three data sources. 
# For tissues with multiple sub-tissues (brain regions, blood cells, lymphoid tissues 
#                                        and intestine) 
# the maximum of all sub-tissues is used for the tissue type. 
# The tab-separated file includes Ensembl gene identifier ("Gene"), 
# analysed sample ("Tissue") and normalized expression ("NX"). 
# The data is based on The Human Protein Atlas version 19.1 and Ensembl version 92.38.	`
repair_genes <- read_tsv("~/hesselberthlab/paper_data/data/other/BER_genes.txt") %>% 
        filter(KEGG_BASE_EXCISION_REPAIR != "> Base excision repair" )
repair_genes <- repair_genes$KEGG_BASE_EXCISION_REPAIR
repair_genes <- c(repair_genes, "RNASEH2A", "RNASEH2B", "RNASEH2C")

refdat <- read_tsv("~/Dropbox (Hesselberth Lab)/HesselberthLab/Images/2019-10-31-RIP-talk/data/rna_tissue_consensus.tsv.zip")
unique(refdat$Tissue)
blood <- tribble(~Tissue, ~names,
                 "T-cells", "T",
                 "NK-cells", "NK",
                 "B-cells", "B",
                 "dendritic cells", "DC",
                 "monocytes", "Mono")

ref_blood <- filter(refdat, Tissue %in% blood$Tissue, 
                    `Gene name` %in% repair_genes) %>%
        left_join(blood)

ref_blood %>% filter(`Gene name` == "UNG") %>%
        ggplot(aes(x = names, y = NX, fill = names)) + 
        geom_bar(stat = 'identity') + 
        coord_flip() + 
        scale_fill_manual(values = colors) + 
        theme_cowplot() + 
        theme(legend.position = 'none') + 
        xlab(NULL) + 
        ylab("Normalized UNG mRNA expression")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/UNG expression PA.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)

ref_blood %>% filter(`Gene name` %in% c("RNASEH2A", "RNASEH2B", "RNASEH2C")) %>%
        ggplot(aes(x = names, y = NX, fill = names)) + 
        geom_bar(stat = 'identity') + 
        coord_flip() + 
        scale_fill_manual(values = colors) + 
        theme_cowplot() + 
        theme(legend.position = 'none') + 
        xlab(NULL) + 
        ylab("Normalized mRNA expression") + 
        facet_grid(rows = vars(`Gene name`))


ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/RNASEH2 expression PA.pdf",
       height = 8, width = 4, units = 'in', useDingbats = F)

ref_blood %>% filter(`Gene name` %in% c("MBD4")) %>%
        ggplot(aes(x = names, y = NX, fill = names)) + 
        geom_bar(stat = 'identity') + 
        coord_flip() + 
        scale_fill_manual(values = colors) + 
        theme_cowplot() + 
        theme(legend.position = 'none') + 
        xlab(NULL) + 
        ylab("Normalized MBD4 mRNA expression")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/MBD4 expression PA.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)

ref_blood %>% filter(`Gene name` %in% c("SMUG1")) %>%
        ggplot(aes(x = names, y = NX, fill = names)) + 
        geom_bar(stat = 'identity') + 
        coord_flip() + 
        scale_fill_manual(values = colors) + 
        theme_cowplot() + 
        theme(legend.position = 'none') + 
        xlab(NULL) + 
        ylab("Normalized SMUG1 mRNA expression")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/SMUG1 expression PA.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)

ref_blood %>% filter(`Gene name` %in% c("APEX1")) %>%
        ggplot(aes(x = names, y = NX, fill = names)) + 
        geom_bar(stat = 'identity') + 
        coord_flip() + 
        scale_fill_manual(values = colors) + 
        theme_cowplot() + 
        theme(legend.position = 'none') + 
        xlab(NULL) + 
        ylab("Normalized APE1 mRNA expression")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/Ape1 expression PA.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)


repair_exp <- rownames_to_column(FetchData(pbmc1, 
                                           vars = c('celltype', repair_genes)), 
                                 "cell_id") %>%
        gather(gene, activity, -cell_id, -celltype) %>% 
        filter(celltype != "Platelet")

repair_exp %>% filter(gene == 'UNG') %>%
        activity_plot(., vline = c(1,2,3)) + 
        xlab("UNG expression")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/UNG expression single cell pbmc1.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)

repair_exp %>% filter(gene %in% c("RNASEH2A", "RNASEH2B", "RNASEH2C")) %>%
        activity_plot(., vline = c(2, 4)) + 
        xlab("mRNA expression") + 
        facet_grid(rows = vars(gene))

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/rnaseh2 expression single cell pbmc1.pdf",
       height = 8, width = 4, units = 'in', useDingbats = F)

repair_exp %>% filter(gene == "SMUG1") %>%
        activity_plot(., vline = c(1,2,3)) + 
        xlab("SMUG1 mRNA expression")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/smug1 expression single cell pbmc1.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)

repair_exp %>% filter(gene == "MBD4") %>%
        activity_plot(., vline = c(1,2, 3)) + 
        xlab("MBD4 mRNA expression")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/mbd4 expression single cell pbmc1.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)

repair_exp %>% filter(gene == "APEX1") %>%
        activity_plot(., vline = c(1,2,3)) + 
        xlab("APEX1 mRNA expression")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/ape1 expression single cell pbmc1.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)

## 10x reference data
# expression matrix was downloaded from here:
#https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k

# Make reference seurat object
ref_pbmc <- CreateSeuratObject(counts = Read10X("../data/other/filtered_gene_bc_matrices/GRCh38/"))
ref_pbmc = NormalizeData(object = ref_pbmc, assay = "RNA", method = "LogNormalize")

# Calculating percent mitochondrial reads per cell 
ref_pbmc[["percent.mt"]] <- PercentageFeatureSet(ref_pbmc, pattern = "^MT-")
ref_pbmc <- FindVariableFeatures(ref_pbmc, selection.method = "vst", nfeatures = 5000)

VlnPlot(ref_pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Decide on features and percent.mt cutoffs by looking at VlnPlot above
ref_pbmc <- subset(ref_pbmc, subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 3000 & 
                        percent.mt < 10)
#load cell type reference data
load("../data/other/pbmc_3k_reference.Rdata")

# Find variable features, calculate clusters, UMAP, TSNE
# Assign celltypes from reference single cell dataset

ref_pbmc <- analyze_seurat(ref_pbmc, pbmc)

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

ref_pbmc <- add_simple_celltype(ref_pbmc, df)
save(ref_pbmc, file = "../data/other/pbmc.4k.seurat.Rdata")

## sc plits for repair genes
repair_exp10x <- rownames_to_column(FetchData(ref_pbmc, 
                                           vars = c('celltype', repair_genes)), 
                                 "cell_id") %>%
        gather(gene, activity, -cell_id, -celltype) %>% 
        filter(celltype != "Platelet")

repair_exp10x %>% filter(activity > 0) %>%
        group_by(gene) %>%
        summarise(n()) -> n_cells

repair_exp10x %>% filter(gene == 'UNG') %>%
        activity_plot(., vline = c(1,2,3)) + 
        xlab("UNG expression")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/UNG expression single cell 4kpbmc10x.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)

repair_exp10x %>% filter(gene %in% c("RNASEH2A", "RNASEH2B", "RNASEH2C")) %>%
        activity_plot(., vline = c(2, 4)) + 
        xlab("mRNA expression") + 
        facet_grid(rows = vars(gene))

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/rnaseh2 expression single cell 4kpbmc10x.pdf",
       height = 8, width = 4, units = 'in', useDingbats = F)

repair_exp10x %>% filter(gene == "SMUG1") %>%
        activity_plot(., vline = c(1,2,3)) + 
        xlab("SMUG1 mRNA expression")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/smug1 expression single cell 4kpbmc10x.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)

repair_exp10x %>% filter(gene == "MBD4") %>%
        activity_plot(., vline = c(1,2, 3)) + 
        xlab("MBD4 mRNA expression")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/mbd4 expression single cell 4kpbmc10x.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)

repair_exp10x %>% filter(gene == "APEX1") %>%
        activity_plot(., vline = c(1,2,3)) + 
        xlab("APEX1 mRNA expression")

ggsave("~/hesselberthlab/projects/10x_haircut/paper_figures/plots/revisions/ape1 expression single cell 4kpbmc10x.pdf",
       height = 4, width = 4, units = 'in', useDingbats = F)








