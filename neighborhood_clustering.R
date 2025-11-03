# neighborhood_clustering.R
# Author: Erin McCaffrey 
# Last updated: 250812
# Overview: This script takes the neighborhood matrices generated on the Dang 
# lab Xenium data. It performs k-means clustering to define cellular neighborhoods. 
# Once the optimal neighborhood configuration is determined, it generates an 
# output with each cell assigned to a particular neighborhood for downstream 
# analysis. 

library(dplyr)
library(tidyr)
library(flowCore)           
library(ggplot2)
library(gplots)
library(RColorBrewer)

##..Step 1: Read in data..##

setwd('/Volumes/T7 Shield/Dang_Lab_Xenium/')
neighbor_data <- read.csv('dang_xenium_neighbor-mat_50px.csv')

##..Step 1.5: Optionally subset the data..##
# subset_cells <- c("Other/distal","Oasl2+","Il1rn+","Gpnmb+","Arg1+","Socs3+",
#                   "Mxd1+","Other","Cd24a+","Cd274+","Gpb2+","Acod1+","Cxcl10+",
#                   "Cd14+","Marcksl1+" )
# neighbor_data<- neighbor_data[neighbor_data$index_subset %in% subset_cells, ]

##..Step 2: Create frequency-normalized version of the matrix..##

# replace any NA with 0s in the counts data
counts_cols <- colnames(neighbor_data[,8:17])
neighbor_data <- neighbor_data %>% 
  mutate(across(all_of(counts_cols), ~ replace_na(.x, 0)))

# normalize to neigborhood total for each cell
neighbor_freqs <- neighbor_data
neighbor_freqs[,8:17] <- neighbor_freqs[,8:17]/neighbor_freqs$total_neighbors

# remove NA values
neighbor_freqs <- neighbor_freqs %>% 
  mutate(across(all_of(counts_cols), ~ replace_na(.x, 0)))

# perform a quick check on the normalization
row_sums <- row_sums <- round(rowSums(neighbor_freqs[, counts_cols], na.rm = TRUE), 12)
table(row_sums_equal_1 <- row_sums == 1)

##..Step 3: Perform clustering and determine value for k..##

# k-means cluster
set.seed(123)
k = 10
neighbor_clusters <- kmeans(neighbor_freqs[, counts_cols], k)

# append cluster annotation to original data
neighbor_freqs$neighborhood_cluster <- neighbor_clusters$cluster

# go through all clusters and get mean frequency of each cell type
cell_types <- counts_cols
hm_allclusters <- matrix(, nrow = k, ncol = length(counts_cols))
for(i in 1:k) {
  temp_mat <- neighbor_freqs[neighbor_freqs[,"neighborhood_cluster"] == i, cell_types]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}
hm_allclusters[is.na(hm_allclusters)] <- 0

# add names to rows and cols
rownames(hm_allclusters) <- paste("neighborhood", 1:k, sep = "")
colnames(hm_allclusters) <- counts_cols 
hm_allclusters

# plot heatmap of all clusters
heatmap.2(hm_allclusters, 
          Colv = T, Rowv = T,
          hclustfun = hclust,
          scale = "column",
          dendrogram = c("both"),
          trace = "none",
          col = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
          density.info = 'none')

##..Step 4: Export the neighborhood cluster-annotated data3..##

write.csv(neighbor_freqs, 'dang_xenium_neighborhoods_r50_k10.csv', row.names = F)

