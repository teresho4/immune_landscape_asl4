library(immunarch)
library(readr)
library(tidyr)
library(DescTools)
library(dplyr)

### IMMUNARCH PREPARATION ###

# bcr - is a concatenated BCR table 'filtered_contig_annotations.csv' for all donors and samples
bcr = read.csv('BCR/bcr_filtered_contig_annotations.csv')
out_dir <- "BCR/immunarch/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cols_10x <- c(
  "barcode", "is_cell", "contig_id", "high_confidence", "length",
  "chain", "v_gene", "d_gene", "j_gene", "c_gene",
  "full_length", "productive",
  "fwr1", "fwr1_nt", "cdr1", "cdr1_nt",
  "fwr2", "fwr2_nt", "cdr2", "cdr2_nt",
  "fwr3", "fwr3_nt", "cdr3", "cdr3_nt",
  "fwr4", "fwr4_nt",
  "reads", "umis",
  "raw_clonotype_id", "raw_consensus_id",
  "exact_subclonotype_id"
)

bcr_clean <- bcr %>%
  group_by(barcode_donor) %>%  filter(all(productive == 'true'), sum(chain == "IGH") == 1, sum(chain %in% c("IGK", "IGL")) == 1) %>%
  ungroup() %>% select(any_of(c("Donor", cols_10x)))

donors <- unique(bcr_clean$Donor)

for (d in donors) {
  df <- bcr_clean[bcr_clean$Donor == d, cols_10x]
  df <- df[order(df$chain, decreasing = FALSE), ]
  write.csv(df, file = file.path(out_dir, paste0(d, "_filtered_contig_annotations.csv")), row.names = FALSE, quote = FALSE, na = "")
}

metadata <- bcr %>% distinct(Sample = Donor, Group, Disease, Annotation, Batch, Tissue) %>%
  mutate(Sample = paste0(Sample, "_filtered_contig_annotations"))
write.table(metadata, file = file.path(out_dir, "metadata.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

### IMMUNARCH UPLOAD ###

immdata <- repLoad(out_dir)
distBCR <- seqDist(immdata$data, .col = 'CDR3.nt', .methods = "hamming")
clustBCR <- seqCluster(immdata$data, distBCR, .perc_similarity = 0.9)
clustBCR_df <- bind_rows(
  lapply(names(clustBCR), function(nm) {
    df <- clustBCR[[nm]]
    donor <- sub("_filtered_contig_annotations$", "", nm)
    df %>% mutate(Donor = donor, Clone = paste0("Clone_", Donor, "_", row_number())) %>%
      select(Barcode, Donor, Clone)
  })
)

meta = bcr[, c('barcode','Donor', 'Disease', 'Annotation')]
meta$Barcode = meta$barcode
meta$barcode = NULL
meta = meta[!duplicated(meta),]
clustBCR_df = clustBCR_df %>% separate_longer_delim(Barcode, delim = ";")
clustBCR_df = merge(clustBCR_df, meta, by=c('Barcode', 'Donor'), all=F)

gini_bcr <- clustBCR_df %>% group_by(Donor, Disease, Annotation, Clone) %>%
  summarise(size = n(), .groups = "drop") %>% group_by(Donor, Annotation, Disease) %>%
  summarise(Gini = Gini(size), .groups = "drop")

