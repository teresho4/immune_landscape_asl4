library(dplyr)
library(MASS)

# data is Seurat object
data = readRDS('PBMC/cd8_t.rds')
df = data@meta.data
df = merge(df, data@reductions$umap@cell.embeddings, by=0)
df$Disease = factor(df$Disease , levels = c('Control', 'ALS4'))
df$Group = factor(df$Group , levels = c('E', 'M', 'L'))

# Two-dimensional kernel density estimation
df <- df %>% group_by(Disease, Group) %>%
  group_modify(~{
    dens <- kde2d(.x$UMAP_1, .x$UMAP_2, n = 200)
    ix <- findInterval(.x$UMAP_1, dens$x)
    iy <- findInterval(.x$UMAP_2, dens$y)
    .x$density <- dens$z[cbind(ix, iy)]
    .x
  }) %>%
  ungroup()

ggplot(df, aes(UMAP_1, UMAP_2, color = density)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(~Disease + Group) +
  theme_classic()
