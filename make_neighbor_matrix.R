# make_neighbor_matrix.R
# Author: Erin McCaffrey 
# Last updated: 250812
# Overview: This script reads in the spatial features for the Dang lab Xenium 
# dataset. For each sample it generates a neighbor matrix with the count of
# neighbors within a selected distance. It creates a version normalized by the 
# total number of neighbors. 

library(dplyr)
library(FNN)
library(spatstat)

##..Step 1: Load the spatial feature data..##

setwd('/Volumes/T7 Shield/Dang_Lab_Xenium/')
spatial_data <- read.csv('Dang-lab_xenium_spatial-data.csv')

##..Step 2: Define parameters..##

radius <- 50
per_sample_results <- list()

##..Step 3: Perform neighbor analysis..##

# Loop over each sample
for (sample_id in unique(spatial_data$Samples)) {
  
  # Subset for the sample
  df_sample <- spatial_data %>% filter(Samples == sample_id)
  
  # Create a spatstat point pattern
  pp <- ppp(
    x = df_sample$x,
    y = df_sample$y,
    window = owin(
      xrange = range(df_sample$x),
      yrange = range(df_sample$y)
    )
  )
  
  # Find all pairs of points within the radius (undirected edges)
  close_pairs <- closepairs(pp, rmax = radius)
  
  # Initialize a list to hold neighbor indices for each point
  neighbor_list <- vector("list", npoints(pp))
  
  # closepairs returns pairs i,j, so fill neighbor_list accordingly
  for (k in seq_along(close_pairs$i)) {
    i <- close_pairs$i[k]
    j <- close_pairs$j[k]
    
    # Add j as neighbor of i
    neighbor_list[[i]] <- c(neighbor_list[[i]], j)
    # Add i as neighbor of j (undirected)
    neighbor_list[[j]] <- c(neighbor_list[[j]], i)
  }
  
  # Get unique cell types in this sample for columns
  cell_types <- unique(df_sample$subset)
  
  # Summarize neighbors' cell types for each cell
  neighbor_summary_counts <- lapply(neighbor_list, function(neigh_idx) {
    if (length(neigh_idx) == 0) {
      setNames(rep(0, length(cell_types)), cell_types)
    } else {
      table(factor(df_sample$subset[neigh_idx], levels = cell_types))
    }
  })
  
  neighbor_counts_df <- do.call(rbind, neighbor_summary_counts)
  total_neighbors <- rowSums(neighbor_counts_df)
  
  # Combine with cell and sample info
  neighbor_summary <- data.frame(
    cell = df_sample$cell,
    Samples = sample_id,
    group = df_sample$Groups,
    index_subset = df_sample$subset,
    index_x = df_sample$x,
    index_y = df_sample$y,
    total_neighbors = total_neighbors,
    neighbor_counts_df,
    row.names = NULL
  )
  
  per_sample_results[[sample_id]] <- neighbor_summary
}

# Combine all sample results
final_neighbor_summary <- do.call(rbind, per_sample_results)

# Export
write.csv(final_neighbor_summary, "dang_xenium_neighbor-mat_50px.csv", row.names = F)
