library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)
library(DescTools)
library(Seurat)
library(purrr)
library(ggpubr)
library(ggalluvial)
`%ni%` = Negate(`%in%`)

############################################## Define clonotype per cell ############################################## 

# Load TCR tables with merged filtered_contig_annotations.csv tables for all samples
# TCR tables are prefilter to cells with annotation from GEX, with exactly one alpha and one beta chains (high-confidence, productive, full-length)
tcr = read.csv('TCR/PBMC_filtered_contig_annotations.csv') # Same for CSF or sorted CD8 T cells
tcr_clones <- tcr %>% mutate(chain_feature = paste(chain, v_gene, j_gene, cdr3_nt, sep = "_")) %>%
  group_by(Donor, barcode_donor, chain) %>%
  summarise(chain_feature = paste(sort(unique(chain_feature)), collapse = "|"),.groups = "drop") %>%
  pivot_wider(names_from = chain, values_from = chain_feature)  %>%
  mutate(clonotype_nt = paste(Donor, TRA, TRB, sep = "_"))

cell_meta <- tcr %>% distinct(barcode_donor, Donor, Group, Disease, Annotation, Major, Batch, Tissue) %>%
  left_join(tcr_clones, by = c("Donor", "barcode_donor")) 

# cell_meta <- cell_meta %>% filter(Major == 'CD8 T')
# write.csv(cell_meta, 'CD8_T_cells_clones.csv')

############################################## Gini ###################################################################

calc_gini <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  if (sum(x) == 0) return(NA_real_)
  Gini(x)
}

# CD4 T or CD8 T
gini_annotation <- cell_meta %>% filter(Major == "CD8 T") %>% 
  count(Donor, Group, Disease, Annotation, clonotype_nt, name = "clone_size") %>%
  group_by(Donor, Group, Disease, Annotation) %>%
  summarise(n_cells = sum(clone_size), n_clones = n(), gini = calc_gini(clone_size), .groups = "drop")

############################################## UMAP ###################################################################

seu = readRDS('PBMC/cd8_t.rds')
seu$barcode_donor = gsub('_.+', '', colnames(seu))
seu$barcode_donor = paste(seu$barcode_donor, seu$Donor, sep='_')

umap_df <- Embeddings(seu, "umap") %>% as.data.frame()
umap_df$barcode_donor = seu$barcode_donor

clone_sizes <- cell_meta %>% filter(Major == "CD8 T") %>%
  count(Donor, clonotype_nt, name = "clone_size")

plot_df <- cell_meta %>% filter(Major == "CD8 T") %>%
  left_join(clone_sizes, by = c("Donor", "clonotype_nt")) %>%
  left_join(umap_df, by = "barcode_donor") 
plot_df$Disease = factor(plot_df$Disease, levels = c('Control', 'ALS4'))
plot_df$Group = factor(plot_df$Group, levels = c('E', 'M', 'L'))

ggplot(plot_df,aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = clone_size), size = 0.1) +
  scale_color_gradient(low = "grey", high = "red") +
  facet_wrap(~Disease+Group) +
  theme_classic() +
  ggtitle("CD8 T cells: clone size")

###################################################### Occupied repertoire space ######################################################

clone_size_df <- cell_meta %>%
  count(Donor, Group, Disease, Major, clonotype_nt, name = "clone_size") %>%
  mutate(
    clone_category = case_when(
      clone_size > 100 ~ "100 < Hyperexpanded",
      clone_size > 20  ~ "20 < Large <= 100",
      clone_size > 5   ~ "5 < Medium <= 20",
      TRUE             ~ "Small <= 5"
    ),
    clone_category = factor(
      clone_category,
      levels = c(
        "Small <= 5",
        "5 < Medium <= 20",
        "20 < Large <= 100",
        "100 < Hyperexpanded"
      )
    )
  )

clone_bar_df <- clone_size_df %>%
  group_by(Donor, Group, Disease, Major, clone_category) %>%
  summarise(cells_in_category = sum(clone_size), .groups = "drop") %>%
  complete(Donor, Group, Disease, Major, clone_category, fill = list(cells_in_category = 0)) %>%
  group_by(Donor, Group, Disease, Major) %>%
  mutate(percent = 100 * cells_in_category / sum(cells_in_category)) %>%
  ungroup()

clone_bar_mean <- clone_bar_df %>%
  group_by(Group, Disease, Major, clone_category) %>%
  summarise(mean_percent = mean(percent, na.rm = TRUE),.groups = "drop")

ggplot(clone_bar_mean, aes(x = Group, y = mean_percent, fill = clone_category)) +
  geom_col(color = "black", width = 0.6) +
  facet_wrap(~Disease) +
  theme_classic() +
  labs(title = "CD8 T cells", x = NULL, y = "Occupied repertoire space (%)", fill = "Clone category")

############################################## TCR sharing for terminal effector ############################################## 

cell_meta =  read.csv('CD4_T_cells_clones.csv')
clones <- cell_meta %>% distinct(Disease, Annotation, clonotype_nt)

overlap_other <- clones %>% rename(Annotation_from = Annotation) %>%
  inner_join(clones %>% rename(Annotation_to = Annotation) %>% dplyr::select(Annotation_to, clonotype_nt),
    by = c("clonotype_nt")) %>% filter(Annotation_to != Annotation_from)

overlap_inside <- clones %>% filter(clonotype_nt %ni% overlap_other$clonotype_nt) %>%
  transmute(Disease, Annotation_from = Annotation, Annotation_to = Annotation, clonotype_nt)

overlap = rbind(overlap_other, overlap_inside) %>% group_by(Disease, Annotation_from, Annotation_to) %>%
  summarise(shared_clones = n_distinct(clonotype_nt), .groups = "drop")
overlap_percent <- overlap %>% group_by(Disease, Annotation_from) %>%
  summarise(Annotation_to, total_clones_from = sum(shared_clones), shared_clones) %>%
  mutate(percent = 100 * shared_clones / total_clones_from) 
overlap_percent$percent = round(overlap_percent$percent, 1)
overlap_percent = overlap_percent %>% filter(percent > 0)
overlap_percent =  overlap_percent %>% filter(Annotation_from == 'Terminal effector')

########################################## CSF clonotypes in PBMC ##########################################

tcr_pbmc = read.csv('TCR/PBMC_filtered_contig_annotations.csv')
tcr_csf = read.csv('TCR/CSF_filtered_contig_annotations.csv')

tcr <- bind_rows(tcr_pbmc, tcr_csf) %>%
  mutate(clonotype = paste(cdr3_nt, v_gene, j_gene, sep = "_"))

## Overall

clones <- tcr %>% filter(Major %in% c("CD4 T", "CD8 T")) %>%
  distinct(Donor, Disease, Major, Tissue, clonotype)

overlap_csf_to_pbmc <- clones %>%
  group_by(Donor, Disease, Major) %>%
  summarise(csf  = list(clonotype[Tissue == "CSF"]), pbmc = list(clonotype[Tissue == "PBMC"]),.groups = "drop") %>%
  rowwise() %>% mutate(n_csf = length(csf), n_shared = length(intersect(csf, pbmc)), 
                       percent_shared = ifelse(n_csf > 0, 100 * n_shared / n_csf, NA_real_)) %>% ungroup()

## By CSF populations

csf_clones <- tcr %>% filter(Major == "CD8 T", Tissue == "CSF") %>%
  distinct(Donor, Disease, Major, Annotation, clonotype)

pbmc_clones <- tcr %>% filter(Major == "CD8 T", Tissue == "PBMC") %>%
  distinct(Donor, Disease, Major, clonotype)

overlap_csf_to_pbmc_annot <- csf_clones %>%
  group_by(Donor, Disease, Major, Annotation) %>%
  summarise(csf = list(unique(clonotype)), .groups = "drop") %>%
  left_join(pbmc_clones %>% group_by(Donor, Disease, Major) %>%
      summarise(pbmc = list(unique(clonotype)), .groups = "drop"), by = c("Donor", "Disease", "Major")) %>%
  rowwise() %>%
  mutate(pbmc = list(if (is.null(pbmc)) character(0) else pbmc), n_csf = length(csf),
    n_shared = length(intersect(csf, pbmc)), percent_shared = ifelse(n_csf > 0, 100 * n_shared / n_csf, NA_real_)) %>% ungroup()

########################################## Cell type correspondece for shared CSF-PBMC clones ##########################################

tcr_pbmc <- read.csv("TCR/PBMC_filtered_contig_annotations.csv")
tcr_csf  <- read.csv("TCR/CSF_filtered_contig_annotations.csv")

tcr_pbmc[tcr_pbmc$Annotation %ni% c("Temra", "HLA-DR+", "Tem GZMB+", "Tem GZMK+"), "Annotation"] = "Other"

tcr <- bind_rows(tcr_pbmc, tcr_csf) %>%
  filter(Major == "CD8 T", chain %in% c("TRA", "TRB")) %>%
  mutate(chain_clonotype = paste(chain, cdr3_nt, v_gene, j_gene, sep = "_"))

tcr <- tcr %>% group_by(Donor, Disease, Tissue, Annotation, barcode) %>%
  summarise(TRA = paste(sort(unique(chain_clonotype[chain == "TRA"])), collapse = "|"),
    TRB = paste(sort(unique(chain_clonotype[chain == "TRB"])), collapse = "|"),
    .groups = "drop") %>% mutate(clonotype = paste(TRA, TRB, sep = "__"))

csf_clones <- tcr %>% filter(Tissue == "CSF")
csf_clones$CSF = csf_clones$Annotation

pbmc_clones <- tcr %>% filter(Tissue == "PBMC")
pbmc_clones$PBMC = pbmc_clones$Annotation

flow_df <- inner_join(csf_clones, pbmc_clones, by = c("Donor", "Disease", "clonotype"))
flow_plot <- flow_df %>% count(Disease, CSF, PBMC, name = "n") %>% group_by(Disease) %>%
  mutate(percent = 100 * n / sum(n)) %>% ungroup()

ggplot(flow_plot, aes(axis1 = CSF, axis2 = PBMC, y = percent )) +
  geom_alluvium(aes(fill = CSF), width = 0.18, alpha = 0.7) +
  geom_stratum(width = 0.18,fill = "white",color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) + facet_wrap(~Disease) + theme_bw() +
  scale_x_discrete(limits = c("CSF", "PBMC"), expand = c(0.1, 0.05) + 
  theme(axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black")) + ylab('% of shared clone pairs')
  
############################################## Public clones ###################################################################

tcr = read.csv('TCR/CSF_filtered_contig_annotations.csv') # CSF, Blood, Sorted CD8 T cells
tcr = tcr %>% filter(chain == 'TRB') # or both chains
tcr_clones <- tcr %>%
  mutate(chain_feature = paste(chain, v_gene, j_gene, cdr3, sep = "_")) %>%
  group_by(Donor, Disease, barcode_donor, chain ) %>%
  summarise(chain_feature = paste(sort(unique(chain_feature)), collapse = "|"),.groups = "drop") %>%
  pivot_wider(names_from = chain, values_from = chain_feature)  %>%
  mutate(clonotype = paste( TRB,  sep = "_"))

cdr3_sharing <- tcr_clones %>% group_by(clonotype) %>%
  summarise(n_donors = n_distinct(Donor), donors = paste(unique(Donor), collapse = ";"),
    donor_disease_status = paste(unique(Disease), collapse = ";"),.groups = "drop")
cdr3_sharing = cdr3_sharing %>% filter(n_donors > 1)

############################################## % overlap VDJdb ###################################################################

vdjdb <- read.csv("VDJdb_HS.tsv", sep = "\t") # Download from https://vdjdb.com/

tcr <- bind_rows(read.csv("TCR/PBMC_filtered_contig_annotations.csv"), read.csv("TCR/SORTED_CD8_filtered_contig_annotations.csv")) %>%
  filter(chain == "TRB") %>% distinct(Donor, Disease, cdr3)

vdjdb_beta <- vdjdb %>% filter(!is.na(CDR3)) %>% distinct(CDR3, Epitope.species, Epitope.gene, Epitope)
total_clones <- tcr %>% group_by(Donor, Disease) %>% summarise(total_beta_clones = n_distinct(cdr3), .groups = "drop")

vdj_hits <- tcr %>% inner_join(bdj_beta, by = c("cdr3" = "CDR3"))
summary <- vdj_hits %>% group_by(Donor, Disease, Epitope.species) %>%
  summarise(vdjdb_hit_clones = n_distinct(cdr3),.groups = "drop") %>%
  left_join(total_clones, by = c("Donor", "Disease")) %>%
  mutate(percent_of_total_beta_clones = 100 * vdjdb_hit_clones / total_beta_clones) %>%
  arrange(Disease, Donor, desc(vdjdb_hit_clones))

