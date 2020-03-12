# Supplementary Fig 3 - qPCR results 

library(tidyverse)
source("scripts/functions.R")

# Load table of delta delta ct for R values
dat <- read_csv("../data/qpcr/deltadeltact_forR.csv")

# Get average and SD from replicate and primers
dat %>% gather(primer, ddct, -Cellline) %>%
        separate(Cellline, into = c("cellline", "rep")) %>%
        filter(primer != "Actin") %>%
        separate(primer, sep = -1, into = c("target", "number")) %>%
        group_by(cellline, target) %>%
        summarize(avg_ddct = mean(ddct),
                  stdev = sd(ddct)) -> plot_dat


# Plot
ggplot(plot_dat, aes(x = cellline, y = avg_ddct, fill = cellline)) + 
        geom_bar(stat = "identity") + 
        cowplot::theme_cowplot() + 
        theme(legend.position = 'none') + 
        xlab(NULL) + 
        ylab("Fold change from Hap1") + 
        scale_fill_manual(values = colors) + 
        geom_errorbar(aes(ymin = avg_ddct-stdev, ymax = avg_ddct+stdev),
                      width = 0.2, position = position_dodge(0.9)) + 
        facet_grid(~target)


#ggsave("../plots/qprc/201908_qpcr_RNASEH2C_UNG.pdf",
#       height =4, width = 8, units = 'in', useDingbats = F)

# Get average and SD from replicate but keep primers separate
dat %>% gather(primer, ddct, -Cellline) %>%
        separate(Cellline, into = c("cellline", "rep")) %>%
        filter(primer != "Actin") %>%
        separate(primer, sep = -1, into = c("target", "number")) %>%
        group_by(cellline, target, number) %>%
        summarize(avg_ddct = mean(ddct),
                  stdev = sd(ddct)) -> plot_dat

# Plot 
ggplot(plot_dat, aes(x = cellline, y = avg_ddct, fill = cellline)) + 
        geom_bar(stat = "identity") + 
        cowplot::theme_cowplot() + 
        theme(legend.position = 'none') + 
        xlab(NULL) + 
        ylab("Fold change from Hap1") + 
        scale_fill_manual(values = colors) + 
        geom_errorbar(aes(ymin = avg_ddct-stdev, ymax = avg_ddct+stdev),
                      width = 0.2, position = position_dodge(0.9)) + 
        facet_grid(target~number)

#ggsave("../plots/qprc/201908_qpcr_RNASEH2C_UNG_primer_set_sep.pdf",
#       height = 6, width = 8, units = 'in', useDingbats = F)

