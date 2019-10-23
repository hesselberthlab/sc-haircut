# Supplementary Figure 5 - QC barplots

## qc
# Does adding hairpins reduce mRNA capture
# 10x data from https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k
# _L001_ fastq was truncated to the first 36 million entries to be comparable to the Dec 2018 pbmc sample
# the full _L001_fastq was use to compare to the Feb 2019 PBMC samples

#V3 data was from https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_v3
# The _L001_ fastq was truncated to the first 35 million reads to be compared to the D21 samples from April 2019

# I ran cellranger 3.0 for all of the data.
library(tidyverse)
library(cowplot)
library(scrunchy)
colors <- c(
        palette_OkabeIto_black,
        scales::brewer_pal(palette = "Paired")(12),
        scales::brewer_pal(palette = "Set1")(9),
        scales::brewer_pal(palette = "Set2")(8),
        scales::brewer_pal(palette = "Dark2")(8)
)
#colors = RColorBrewer::brewer.pal(9, "Greys")[c(3, 5, 7, 9, 2)]

sample = list.files(path ="../data/qc/csv", full.names = T)

dat = tibble(sample) %>%
        mutate(conetents = map(sample, read_csv)) %>%
        mutate(sample = str_remove(string = sample, 
                                   pattern = "data/qc/csv/"),
               sample = str_remove(sample, "_metrics_summary.csv"))

dat = unnest(dat)

qc_bar = function(df, var, xlabs  = c("PBMC", "PBMC + substrate")){
        p <- ggplot(df, aes_string(x = "sample", y = var, fill = "sample")) + 
                geom_bar(stat = "identity") + 
                scale_fill_manual(values = colors) + 
                theme_cowplot() + 
                theme(legend.position = "none") + 
                scale_x_discrete(labels= xlabs) + 
                labs(y = str_replace_all(var, "_", " "),
                     x = '')
        p
        
}

colnames(dat) <- str_replace_all(colnames(dat), " ", "_")
dat <- mutate(dat, Sequencing_Saturation = str_remove(Sequencing_Saturation, "%"),
              Sequencing_Saturation = as.numeric(Sequencing_Saturation))

pbmc1 <- filter(dat, sample %in% c("pbmc_36M_from10x", "pbmc1"))

plot_var <- colnames(dat)[c(5, 2, 3, 19, 4, 20, 7)]


# Plotting PBMC1 (Figure 2) QC
plot_df <- tibble(plot_var) %>%
        mutate(plot = map(plot_var, ~ qc_bar(pbmc1, .x, c("PBMC", "70 nM"))))

plot_grid(plotlist = plot_df$plot, nrow = 1)

#ggsave("plots/qc/pbmc1_and_36M_from_10x.pdf",
#       height = 3, width = 15, units = 'in', useDingbats = F)

# Plot supplementary fig 4 qc with 10 nM substrates
pbmc3 <- filter(dat, sample %in% c("pbmc_190M_from10x", "pbmc2", "pbmc3"))
plot_df <- tibble(plot_var) %>%
        mutate(plot = map(plot_var, ~ qc_bar(pbmc3, .x, c("PBMC", "10 nM ",
                                                          "10 nM "))))

plot_grid(plotlist = plot_df$plot, nrow = 1)
#ggsave("plots/qc/pbmc2and3_and_190M_from_10x.pdf",
#       height = 3, width = 15, units = 'in', useDingbats = F)


## Time course/dilution experiment
# QC 
# Reference data from https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v2
# PBMC 60 sample was downsampled to 10 million reads to match the sequencing depth of PBMC15 and 30

pbmc <- filter(dat, files %in% c("10x_pbmc", "pbmc_15", "pbmc_30", "pbmc_60_10mllion"))

plot_var <- colnames(dat)[c(5, 2, 3, 19, 4, 20, 7)]

plot_df <- tibble(plot_var) %>%
        mutate(plot = map(plot_var, ~ qc_bar(pbmc, .x, c("10x", "15", "30", "60"))))

plot_grid(plotlist = plot_df$plot, nrow = 1)
#ggsave("plots/pbmc/qc_plots.pdf",
#       height = 3, width = 15, units = "in", useDingbats = F)


