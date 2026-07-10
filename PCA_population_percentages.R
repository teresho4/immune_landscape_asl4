library(dplyr)
library(tidyverse)

# load tables with percentages from each major population
percentages = rbind(percentages_cd8, percentages_cd4, percentages_gt, percentages_b, percentages_nk, percentages_myeloid)
percentages = percentages[, c('Donor', 'Annotation', 'percent')]
df_wide <- percentages %>% pivot_wider(names_from = Annotation, values_from = percent)
df_wide <- as.data.frame(df_wide)
row.names(df_wide) <- df_wide$Donor
df_wide$Donor <- NULL

pca <- prcomp(df_wide, scale = T, center = T)

# Make PCA data frame
pca_df <- as.data.frame(pca$x)
meta <- data@meta.data[, c('Donor', 'Source')]
meta <- meta[!duplicated(meta),]
pca_df$Donor = row.names(pca_df)
pca_df <- left_join(pca_df, meta, by = "Donor")

ggplot(pca_df, aes(x = PC1, y = PC2, color=Source)) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal() 

