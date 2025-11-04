# neighborhood_functional_markers.R
# Author: Erin McCaffrey 
# Last updated: 250930
# Overview: This script takes the neighborhood data from the Dang lab. It 
# appends the neighborhood assignment to the metadata then evaluates expression
# of certain genes across all neighborhoods. 

library(dplyr)
library(ggplot2)
library(gplots)
library(reshape2)
library(RColorBrewer)


##...Load in data..##

setwd('/Volumes/T7 Shield/Dang_Lab_Crypto_Grans/xenium')
neighborhood_data <- read.csv('dang_xenium_neighborhoods_r50_k10.csv')
neighborhood_color_key <- read.csv('r50_k10_color-key.csv')
spatial_data <- read.csv('Dang-lab_xenium_spatial-data.csv')

##..Add the neighborhood assignment to the spatial data..##

# subset neighborhood data
neighborhood_data_sub <- neighborhood_data %>%
  select(cell, index_subset, neighborhood_cluster)

# merge with spatial data
spatial_data_neighborhood <- merge(spatial_data, neighborhood_data_sub, by = c('cell'))

##..Summarize mean expression per gene per neighborhood..##

# define genes
genes_of_interest <- c(
  "Mxd1",
  "Cd274",
  "Arg1",
  "Chil3",
  "Cd68",
  "Itgam",
  "Ighm",
  "Tcf7",
  "Cd3e",
  "Col1a2",
  "Cxcl12",
  "Cxcl10",
  "Cxcl16",
  "Cxcl13")

# summarize
by_cluster_gene_means <- spatial_data_neighborhood %>%
  group_by(neighborhood_cluster) %>%
  summarise(across(all_of(genes_of_interest), ~ mean(.x, na.rm = TRUE)),.groups = "drop")

##..Generate a heatmap..##

# define order of clusters to match cell subset heatmap
cluster_order <- c(4,7,6,2,9,5,1,10,3,8)
hmap_data <- by_cluster_gene_means %>% 
  mutate(neighborhood_cluster = factor(neighborhood_cluster,
                                     levels = c(4,7,6,2,9,5,1,10,3,8),
                                     ordered = TRUE)) %>%
  arrange(neighborhood_cluster)

# create matrix
hmap_data.m <- as.matrix(hmap_data[,-1])
rownames(hmap_data.m) <- hmap_data$neighborhood_cluster

# plot
heatmap.2(hmap_data.m, 
          Colv = T, Rowv = F,
          hclustfun = hclust,
          scale = "column",
          dendrogram = c("column"),
          trace = "none",
          col = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
          density.info = 'none')

##..Plot violin of expression across clusters..##

# melt the data 
data.m <- melt(spatial_data_neighborhood, id.vars = c('Samples','Groups','index_subset','neighborhood_cluster'),
               measure.vars = genes_of_interest)

# get neighborhood colors for plotting
neighborhood_clusters <-  levels(factor(data.m$neighborhood_cluster))
plot_colors <- droplevels(neighborhood_color_key[neighborhood_color_key$neighborhood_cluster %in% neighborhood_clusters,])
plot_colors$neighborhood_cluster <- factor(plot_colors$neighborhood_cluster, levels = neighborhood_clusters)
plot_colors <- plot_colors[order(plot_colors$neighborhood_cluster),]
color <- as.vector(plot_colors$hex)

# create violins (1 per marker)

ggplot(data.m, aes(x = neighborhood_cluster, y = value, fill =as.factor(neighborhood_cluster))) +
  geom_violin(trim = TRUE) + 
  scale_fill_manual(values = color) +
  facet_wrap(~variable, scales = "free_y", ncol = 7) 
  

