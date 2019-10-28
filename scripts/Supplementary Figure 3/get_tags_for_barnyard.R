## Pull out cell barcodes for each cell type for mRNA coverage plots

load("../data/barnyard/barnyard.seurat.object.Rdata")

#Getting hairpin info from seurat object

df <- rownames_to_column(as.data.frame(t(as.matrix(GetAssayData(barnyard_seurat, assay = 'repair', slot = 'counts')))), "cell_id") %>%
        gather(hairpin_pos, count, -cell_id) %>%
        separate(hairpin_pos, into = c("hairpin", "position")) %>%
        mutate(position = as.double(position),
               count = as.double(count))


repair_position = data_frame(hairpin = c('Uracil', 'riboG'),
                             repair_position = c(45, 44))

df %>% left_join(repair_position) %>%
        filter(position == repair_position) -> rt

rt %>% group_by(hairpin, position) %>%
        summarize(max_count = max(count)) %>%
        mutate(cut_off = round(max_count * .05)) -> cutoffs

r = cutoffs$cut_off[1]
u = cutoffs$cut_off[2]

rt %>% select(-position, -repair_position) %>% 
        spread(hairpin, count) %>%
        mutate(color = if_else(riboG >= r & Uracil <+ u, 'UNGKO',
                               if_else(Uracil >= u & riboG <= r, 'RNASEH2KO',
                                       if_else(riboG >= r & Uracil >= u, 'Both',
                                               'Low signal')))) %>% 
        select(-riboG, -Uracil) -> color_df



color_df %>% filter(color == "RNASEH2KO") %>% 
        mutate(tag = paste0("CB:Z:", cell_id, "-1")) %>% 
        select(tag) %>%
        write_tsv(path = "../data/barnyard_102018_RNASEH2CKO_cell_tags.tsv",
                  col_names = F)

color_df %>% filter(color == "UNGKO") %>% 
        mutate(tag = paste0("CB:Z:", cell_id, "-1")) %>% 
        select(tag) %>%
        write_tsv(path = "../data/barnyard_102018_UNGKO_cell_tags.tsv",
                  col_names = F)

color_df %>% filter(color %in% c("UNGKO", "RNASEH2KO")) %>%
        mutate(tag = paste0("CB:Z:", cell_id, "-1")) %>% 
        select(tag, color) %>%
        write_tsv(path = "../data/barnyard_102018_celltype_tags.tsv",
                  col_names = F)
