### First run arcasHLA in bash for each donor-specific BAM file ###
# arcasHLA extract ${bam} -o ${OUTDIR} -t 8 -v
# arcasHLA genotype "hla_output/$sample/${sample}.extracted.fq.gz" -o "hla_output/$sample" -t 8 -v (all genes by default)

library(jsonlite)
library(dplyr)

files = list.files('.')
hla = NULL
for(file in files){
  data <- fromJSON(paste(file, '/', file,'.genotype.json', sep=''))
  df <- stack(data)
  colnames(df) <- c("Allele", "Gene")
  df$Donor = file
  hla = rbind(hla, df)
}

meta = read.csv('metadata_donor.csv', row.names = 'X')
hla = merge(hla, meta, by='Donor')

# Collapse alleles per gene
haplotypes <- hla %>%
  group_by(Donor, Disease, Gene) %>%
  summarise(Alleles = paste(unique(Allele), collapse = " / "), .groups = "drop") %>%
  pivot_wider(names_from = Gene, values_from = Alleles)

# Define known combinations
group_map <- list(
  # DQ
  "DQ2.5" = list(DQA1="DQA1*05:01", DQB1="DQB1*02:01"),
  "DQ2.2" = list(DQA1="DQA1*02:01", DQB1="DQB1*02:02"),
  "DQ8"   = list(DQA1="DQA1*03:01", DQB1="DQB1*03:02"),
  "DQ7"   = list(DQA1="DQA1*05:05", DQB1="DQB1*03:01"),
  "DQ9"   = list(DQA1="DQA1*03:02", DQB1="DQB1*03:03"),
  
  # DR-DQ linkage blocks
  "DR3-DQ2.5" = list(DRB1="DRB1*03:01", DQA1="DQA1*05:01", DQB1="DQB1*02:01"),
  "DR4-DQ8"   = list(DRB1="DRB1*04:01", DQA1="DQA1*03:01", DQB1="DQB1*03:02"),
  "DR7-DQ2.2" = list(DRB1="DRB1*07:01", DQA1="DQA1*02:01", DQB1="DQB1*02:02"),
  
  # DP groups
  "DP2.1" = list(DPA1="DPA1*02:01", DPB1="DPB1*01:01"),
  "DP4"   = list(DPA1="DPA1*01:03", DPB1="DPB1*04:01"),
  "DP5"   = list(DPA1="DPA1*01:03", DPB1="DPB1*02:01")
)

# Match donor haplotypes to known groups
assign_groups <- function(row, group_map) {
  matched <- c()
  
  for (grp in names(group_map)) {
    required <- group_map[[grp]]
    ok <- all(sapply(names(required), function(gene) {
      if (!gene %in% names(row)) return(FALSE)
      val <- as.character(row[[gene]])
      if (is.na(val)) return(FALSE)
      str_detect(val, fixed(required[[gene]]))
    }))
    if (ok) matched <- c(matched, grp)
  }
  if (length(matched) == 0) {
    return(NA_character_)
  }
  paste(matched, collapse = "; ")
}

haplotypes <- haplotypes %>% rowwise() %>% mutate(HLA_Group = assign_groups(cur_data(), group_map)) %>% ungroup()
haplotypes[is.na(haplotypes$HLA_Group), 'HLA_Group'] = 'No matching group'
