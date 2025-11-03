# neighborhood_cluster_visualization.R
# Author: Erin McCaffrey 
# Last updated: 250813
# Overview: This script takes the clustered neighborhood data and visualizes
# the clusters in the spatial positions of the Xenium data.

library(Seurat)
library(dplyr)
library(purrr)
library(ggrastr)
library(ggplot2)
library(tibble)

##..Step 1: Load the RDS file and the neighborhood clusters..##

setwd("/Volumes/T7 Shield/dang_lab_crypto_grans/Xenium/")
neighborhood_data <- read.csv('dang_xenium_neighborhoods_r50_k10.csv')
neighborhood_color_key <- read.csv('r50_k10_color-key.csv')

##..Step 2: Define the output directory and sample list..##
out_dir <- "/Volumes/T7 Shield/dang_lab_crypto_grans/Xenium/plots"
samples <- sort(unique(neighborhood_data$Samples))

##..Step 3: Visualize the neighborhood clusters, plot and save..##

# (one sample at a time)

plot_data <- neighborhood_data[neighborhood_data$Samples == 'd42_2',]
neighborhood_clusters <- unique(plot_data$neighborhood_cluster)
plot_colors <- droplevels(neighborhood_color_key[neighborhood_color_key$neighborhood_cluster %in% neighborhood_clusters,])
plot_colors$Phenotype <- factor(plot_colors$neighborhood_cluster, levels = neighborhood_clusters)
plot_colors<-plot_colors[order(plot_colors$neighborhood_cluster),]
color<-as.vector(plot_colors$hex)

# ggplot(plot_data, aes(x = index_x, y = index_y, color = as.factor(neighborhood_cluster))) +
#   geom_point(size = 0.25) +
#   scale_color_manual(values = color) +
#   coord_fixed() +
#   theme_void() +
#   theme(legend.position = "right") +
#   labs(color = "Cluster")

ggplot(plot_data, aes(x = index_x, y = index_y, color = as.factor(neighborhood_cluster))) +
  geom_point_rast(size = 0.25, raster.dpi = 300) +
  scale_color_manual(values = color) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "right") +
  labs(color = "Cluster")


# (loop over all samples)

walk(samples, function(sid) {
  # subset to one sample
  plot_data <- neighborhood_data %>% filter(Samples == sid)
  if (nrow(plot_data) == 0) return(invisible(NULL))
  
  # clusters present in this sample (keep this order)
  neighborhood_clusters <- unique(plot_data$neighborhood_cluster)
  
  # color key for only the present clusters, ordered to match the vector above
  plot_colors <- neighborhood_color_key %>%
    filter(neighborhood_cluster %in% neighborhood_clusters) %>%
    distinct(neighborhood_cluster, hex) %>%
    mutate(order = match(neighborhood_cluster, neighborhood_clusters)) %>%
    arrange(order)
  
  # named color vector for scale_color_manual
  color_vec <- setNames(plot_colors$hex, plot_colors$neighborhood_cluster)
  
  p <- ggplot(
    plot_data,
    aes(x = index_x, y = index_y, color = as.factor(neighborhood_cluster))
  ) +
    geom_point_rast(size = 0.05, raster.dpi = 300) +
    scale_color_manual(
      values = color_vec,
      breaks = neighborhood_clusters,  # keep legend order consistent
      drop = FALSE
    ) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "right") +
    labs(color = "Cluster", title = sid)
  
  # save (PNG keeps the rasterization; use PDF if you prefer vectors elsewhere)
  outfile <- file.path(out_dir, paste0("neighborhood_", sid, ".png"))
  ggsave(outfile, p, width = 10, height = 10, dpi = 300, bg = "white")
})
