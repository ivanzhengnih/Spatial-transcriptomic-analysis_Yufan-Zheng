# xenium_extract_cell_spatial_features.R
# Author: Erin McCaffrey 
# Date created: 250811
# Overview: This script reads in Ivan's Xenium data. It extracts out the position
# and identity of each cell. Additionally, it defines the granuloma x and y center
# for each tissue

library(Seurat)
library(dplyr)
library(spacexr)

##..Step 1: Load the RDS file..##

setwd('/Volumes/ThunderBlade/Dang_Lab_Xenium')
seurat_object <- readRDS("Xenium_Ivan_merged_clustered.rds")
seurat_macs <- readRDS("Xenium_Ivan_Macrophages.rds")

##..Step 2: Pull out the cell data..##

# extract the metadata
metadata <- seurat_object@meta.data
metadata$cell <- rownames(metadata)

##..Step 3: Pull out the spatial information..##

# extract the spatial coordinates
spatial_coords  <- bind_rows(lapply(Images(seurat_object), function(x) GetTissueCoordinates(seurat_object, image = x, scale = "hires")))
spatial_coords_centroid <- aggregate(. ~ cell, data = spatial_coords, FUN = mean)

# merge the metadata and spatial data
metadata_spatial_coords <- merge(metadata, spatial_coords_centroid, by = c('cell'))

##..Step 4: Append the macrophage subsets..##

# get just the mac cell ids and subset name
metadata_macs <- seurat_macs@meta.data
metadata_macs$cell <- rownames(metadata_macs)
keep_cols <- c('cell','cell_types')
metadata_macs_subset <- metadata_macs %>% select(all_of(keep_cols))

# add a new cluster id column to the metadata
metadata_spatial_coords$subset <- metadata_spatial_coords$major_cell_types_no_spaces
merged_df<- merge(metadata_spatial_coords, metadata_macs_subset, by = "cell", all.x = TRUE)
merged_df <- merged_df %>%
  mutate(subset = if_else(!is.na(cell_types), cell_types, subset))
merged_df <- merged_df %>% select(-cell_types)
merged_df$subset <- droplevels(merged_df$subset)

##..Step 5: Export..##
write.csv(merged_df, file="Dang-lab_xenium_spatial-data.csv",row.names = FALSE)
