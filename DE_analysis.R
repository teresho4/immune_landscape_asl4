library(Seurat)
library(dplyr)
library(purrr)
library(presto)
library(msigdbr)
library(fgsea)
library(edgeR)
library(ggplot2)
library(limma)

########################## DE LOW-EXPANDED VS HIGH-EXPANDED CLONES ##########################

cd8_seurat = readRDS('PBMC/cd8_t.rds')
# table with clone size for each cell in CD8 T cells
tcr_clones <- read.csv('CD8_T_cells_clones.csv') 
tcr_clones <- tcr_clones %>% mutate(ExpansionGroup = case_when(Clone_size <= 5 ~ "LowExpanded", 
                                                               Clone_size > 20 ~ "HighExpanded", TRUE ~ NA_character_))
tcr_clones = tcr_clones[!is.na(tcr_clones$ExpansionGroup),]

cd8_seurat = subset(cd8_seurat, barcode_donor %in% clone_meta$barcode_donor)
cd8_seurat$ExpansionGroup = "LowExpanded"
tcr_clones = tcr_clones %>% filter(ExpansionGroup == "HighExpanded")
cd8_seurat@meta.data[cd8_seurat$barcode_donor %in% tcr_clones$barcode_donor, "ExpansionGroup"] = "HighExpanded"

genes = row.names(cd8_seurat)
genes <- genes[!grepl("^(TRA|TRB|IGH|IGK|IGL|MT-|RPL|RPS)", genes)]

data = subset(cd8_seurat, Annotation == 'Tem GZMB+') # select here the population of interest
result <- list()
for (donor in unique(data$Donor)) {
  obj <- subset(data, Donor == donor)
  de <- tryCatch({
    wilcoxauc(obj, group_by = "ExpansionGroup", seurat_assay = "RNA", assay = "data") %>%
      filter(feature %in% genes, group == "HighExpanded") %>%
      mutate(Donor = donor, contrast = "HighExpanded_vs_LowExpanded")
  }, error = function(e) {
    data.frame()
  })
  result[[donor]] <- de
}

meta_results <- bind_rows(result, .id = "Donor") %>%
  rename(gene = feature) %>%
  filter(!is.na(logFC), !is.na(pval)) %>%
  mutate(z = abs(qnorm(pval / 2)), SE = abs(logFC / z)) %>%
  filter(is.finite(SE), SE > 0) %>%
  group_by(gene) %>%
  summarise(meta = list(
      tryCatch(
        rma.uni(yi = logFC, sei = SE, method = "REML"), 
        error = function(e) NULL)),.groups = "drop")

meta_table <- meta_results %>%
  filter(!sapply(meta, is.null)) %>%
  mutate(
    meta_logFC = sapply(meta, function(x) x$b[1]),
    meta_SE    = sapply(meta, function(x) x$se),
    meta_z     = sapply(meta, function(x) x$zval),
    meta_p     = sapply(meta, function(x) x$pval),
    ci_lower   = sapply(meta, function(x) x$ci.lb),
    ci_upper   = sapply(meta, function(x) x$ci.ub),
    tau2       = sapply(meta, function(x) x$tau2)
  ) %>%
  select(-meta) %>%
  mutate(p_adj = p.adjust(meta_p, method = "BH"))

########################## DE ALS4 VS CONTROL ##########################

data = readRDS('PBMC/pbmc_filtered.rds')
data = subset(data, Source == 'NIH')
data = subset(data, Batch == 'Our_1')

meta <- data@meta.data[, c('Donor', 'Disease')]
meta = meta[!duplicated(meta),]
row.names(meta) = meta$Donor

single_cell_result = list()
for(subpopulation in unique(data$Annotation_subpopulations)){
  sub = subset(data, Annotation_subpopulations == subpopulation)
  counts <- AggregateExpression(sub, assays = "RNA", group.by = "Donor", slot = "counts", return.seurat = FALSE)$RNA
  counts = as.data.frame(counts)
  meta <- data@meta.data[, c('Donor', 'Disease')]
  meta = meta[!duplicated(meta),]
  row.names(meta) = meta$Donor
  meta = meta[colnames(counts),]
  tryCatch(
    {
      dge <- DGEList(counts = counts)
      keep <- filterByExpr(dge, group = meta$Disease)
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      dge <- calcNormFactors(dge)
      
      design <- model.matrix(~ Disease, data = meta)
      v <- voom(dge, design, plot = TRUE)
      fit <- lmFit(v, design)
      fit <- eBayes(fit)
      result <- topTable(fit, coef = "DiseaseControl", number = Inf, sort.by = "P")
    }, error = function(e) {result <- NULL}
  )
  single_cell_result[[subpopulation]] <- result
}

# So that values are relative to ALS4
for(name in names(single_cell_result)){
  single_cell_result[[name]]$logFC  = single_cell_result[[name]]$logFC * -1
  single_cell_result[[name]]$t  = single_cell_result[[name]]$t * -1
}

de_counts <- imap_dfr(single_cell_result, ~{
  tibble(
    Cluster = .y,
    Direction = c("Up", "Down"),
    N = c(
      sum(.x$adj.P.Val < 0.05 & .x$logFC > 0, na.rm = TRUE),
      sum(.x$adj.P.Val < 0.05 & .x$logFC < 0, na.rm = TRUE)
    )
  )
}) 

de_counts[de_counts$Direction == 'Down', 'N'] = de_counts[de_counts$Direction == 'Down', 'N'] * -1
de_counts <- de_counts %>% group_by(Cluster) %>% mutate(total = sum(abs(N))) %>% ungroup()
de_counts <- de_counts %>% filter(total > 5)

ggplot(de_counts, aes(x = reorder(Cluster, total), y = N, fill = Direction)) +
  geom_col(color='black') +
  coord_flip()+ theme_bw() + xlab('') + ylab('Number of DE genes')+ 
  scale_y_continuous(labels = c("-1000" = "1000", "-500" = "500", "0" = "0", "500" = "500","1000" = "1000"))+ 
  theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black")) + 
  scale_fill_manual(values = c("#5d7994", "#c4401b"))

# Volcano for C mono
df = single_cell_result$`Myeloid C mono`
ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = adj.P.Val < 0.05)) +
  ggrastr::geom_point_rast() + theme_bw() + 
  ggrepel::geom_text_repel(data = df %>% filter(-log10(adj.P.Val) > 2.8), aes(label=gene), max.overlaps = 15)+ 
  scale_color_manual(values = c('grey', '#bf5252')) + ggtitle('Classical monocytes, ALS4 vs Control')

## Enrichment plots ##

# Current study to Zhang et al.
gene_set = list()
gene_set[['ALS4']] = single_cell_result$`Myeloid C mono`[order(single_cell_result$`Myeloid C mono`$t, decreasing = T), 'gene'][1:100]
public = read.csv('de_paper.csv', sep='\t') # Signature from Zhang et al.
rank = public$avg_log2FC
names(rank) = public$gene_id
fgseaMultilevel(gene_set, rank)
plotEnrichment(gene_set$ALS4, rank)
# Zhang et al. to current study
rank = single_cell_result$`Myeloid C mono`$t
names(rank) = single_cell_result$`Myeloid C mono`$gene
gene_set = list()
gene_set[['ALS4']] = public[order(public$avg_log2FC, decreasing = T), 'gene_id'][1:100]
fgseaMultilevel(gene_set, rank)
plotEnrichment(gene_set$ALS4, rank)

# HALLMARK
gene_sets = msigdbr(species = "Homo sapiens")
gene_sets  <- gene_sets %>% dplyr::filter(gs_cat == "H")
msigdbr_list = split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
fgseaResults_sample <- fgseaMultilevel(msigdbr_list, rank, minSize = 15, maxSize = 500)
plotEnrichment(msigdbr_list$HALLMARK_TNFA_SIGNALING_VIA_NFKB, rank)
plotEnrichment(msigdbr_list$HALLMARK_INFLAMMATORY_RESPONSE, rank)                                                                                                 
                                                                                                                                                                   
## PCA for pseudobulk expression PBMC-level ##
data = readRDS('PBMC/pbmc_filtered.rds')
counts <- AggregateExpression(data, assays = "RNA", group.by = "Donor", slot = "counts", return.seurat = FALSE)$RNA
dge <- DGEList(counts)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE, prior.count = 1)
pca <- prcomp(t(logCPM), scale. = TRUE)
pca_df <- data.frame(Donor = colnames(logCPM), PC1 = pca$x[,1], PC2 = pca$x[,2])
meta = data@meta.data[, c('Donor', 'Disease', 'Source')]
meta = meta[!duplicated(meta),]
pca_df = merge(pca_df, meta, by='Donor')
ggplot(pca_df, aes(PC1, PC2, color=Disease)) +
  geom_point(size = 3) +
  theme_classic()
