# xenium_color_segmentation_masks.R
# Author: Erin McCaffrey 
# Date created: 250930
# Overview: This script reads in the Seurat object, the annotated metadata, 
# and user-defined color keys for the cell subsets and the neighborhood 
# clusters. It creates images of each tissue with the cells colored by 
# subset and neighborhood cluster assignment.

library(Seurat)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(rlang)

##..Load data file..##

setwd('/Volumes/T7 Shield/Dang_Lab_Crypto_Grans/xenium')
seurat_object <- readRDS("./xenium_data/Xenium_Ivan_merged_clustered.rds")
neighborhood_data <- read.csv('dang_xenium_neighborhoods_r50_k10.csv')
spatial_data <- read.csv('Dang-lab_xenium_spatial-data.csv')
neighborhood_color_key <- read.csv('r50_k10_color-key.csv')
cell_color_key <- read.csv('major_subsets_color-key.csv')
mac_color_key <- read.csv('mac_subset_color-key.csv')

##..Merge the neighborhood data with the metadata..##

# subset neighborhood data
neighborhood_data_sub <- neighborhood_data %>%
  select(cell, index_subset, neighborhood_cluster)

# merge with spatial data
spatial_data_neighborhood <- merge(spatial_data, neighborhood_data_sub, by = c('cell'))

##..Define the plotting attributes..##

# data
ann_df <- spatial_data_neighborhood %>% select(cell, major_cell_types) %>% distinct()
anno_col  <- "major_cell_types" 

# convert color key to named vector 
label_col <- "major_cell_types"
hex_col   <- "hex"

color_key_clean <- cell_color_key %>%
  transmute(
    label = trimws(as.character(.data[[label_col]])),
    hex   = trimws(as.character(.data[[hex_col]]))
  ) %>%
  filter(!is.na(label), label != "", !is.na(hex), hex != "") %>%
  distinct(label, .keep_all = TRUE)   # keep first if duplicates

annotation_colors <- setNames(color_key_clean$hex, color_key_clean$label)

##..Plot..##

# output dir
out_dir <- "xenium_segmentation_cell_pngs"
dir.create(out_dir, showWarnings = FALSE)

# Helper functions: 

# get segmentation for one image and normalize columns 
get_seg_one <- function(obj, im) {
  df <- GetTissueCoordinates(obj, image = im, type = "segmentations")
  df$image <- im
  
  nm <- names(df)
  # Ensure we have cell / x / y column names
  if (!"cell" %in% nm) {
    if ("cell_id" %in% nm) names(df)[match("cell_id", nm)] <- "cell"
    else if ("barcode" %in% nm) names(df)[match("barcode", nm)] <- "cell"
    else stop("No 'cell' column found in segmentation table for image: ", im)
  }
  if (!"x" %in% nm && "vertex_x" %in% nm) names(df)[match("vertex_x", nm)] <- "x"
  if (!"y" %in% nm && "vertex_y" %in% nm) names(df)[match("vertex_y", nm)] <- "y"
  
  # ring/vertex ordering
  if (!"ring_id" %in% names(df))      df$ring_id <- 1L
  if (!"vertex_order" %in% names(df)) df <- df %>% group_by(cell, ring_id) %>% mutate(vertex_order = row_number()) %>% ungroup()
  
  df %>% arrange(cell, ring_id, vertex_order)
}

# get segmentation for one image and normalize columns 
add_scalebar <- function(p, df, length_um = 50, flip_y = FALSE,
                         corner = c("ll","lr","ul","ur"),
                         pad_frac = 0.02, bar_height_frac = 0.006, size_text = 3) {
  corner <- match.arg(corner)
  xr <- range(df$x, na.rm = TRUE); yr <- range(df$y, na.rm = TRUE)
  xpad <- diff(xr) * pad_frac; ypad <- diff(yr) * pad_frac
  bar_h <- diff(yr) * bar_height_frac
  
  # x position
  if (corner %in% c("lr","ur")) { x1 <- xr[2] - xpad; x0 <- x1 - length_um
  } else { x0 <- xr[1] + xpad;   x1 <- x0 + length_um }
  
  # y position (respecting flip_y)
  if ((!flip_y && corner %in% c("ll","lr")) || (flip_y && corner %in% c("ul","ur"))) {
    # bottom
    y0 <- if (!flip_y) yr[1] + ypad else yr[2] - ypad
    y_text <- y0 + if (!flip_y) bar_h * 2.2 else -bar_h * 2.2
    if (flip_y) bar_h <- -bar_h
  } else {
    # top
    y0 <- if (!flip_y) yr[2] - ypad else yr[1] + ypad
    y_text <- y0 - if (!flip_y) bar_h * 2.2 else -bar_h * 2.2
    if (!flip_y) bar_h <- -bar_h
  }
  
  p +
    geom_rect(aes(xmin = x0, xmax = x1, ymin = y0, ymax = y0 + bar_h),
              inherit.aes = FALSE, fill = "black") +
    annotate("text", x = (x0 + x1)/2, y = y_text,
             label = paste0(length_um, " µm"),
             vjust = 0.5, hjust = 0.5, size = size_text)
}

# Iterate images, plot & save 
img_names <- Images(seurat_object)
flip_y <- FALSE  

for (i in seq_along(img_names)) {
  im <- img_names[i]
  message(sprintf("Rendering %s (%d/%d)", im, i, length(img_names)))
  
  seg <- get_seg_one(seurat_object, im)
  
  plot_df <- seg %>%
    left_join(ann_df, by = "cell") %>%
    mutate(!!anno_col := tidyr::replace_na(.data[[anno_col]], "Unlabeled"))
  
  # Ensure color key covers all labels we see
  cols <- annotation_colors
  unseen <- setdiff(unique(plot_df[[anno_col]]), names(cols))
  if (length(unseen)) {
    warning("Annotations not in color key: ", paste(unseen, collapse = ", "),
            " — assigning '#C0C0C0'.")
    cols <- c(cols, setNames(rep("#C0C0C0", length(unseen)), unseen))
  }
  plot_df[[anno_col]] <- factor(plot_df[[anno_col]], levels = names(cols))
  
  p <- ggplot(plot_df, aes(x = x, y = y,
                           group = interaction(cell, ring_id),
                           fill  = .data[[anno_col]])) +
    geom_polygon(color = NA) +            # set color="white", linewidth=0.05 for borders
    coord_equal()
  
  if (flip_y) p <- p + scale_y_reverse()
  
  p <- p +
    scale_fill_manual(values = cols, drop = FALSE) +
    labs(x = NULL, y = NULL, fill = anno_col, title = im) +
    theme_void(base_size = 12) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))
  
  # optional: add scale bar
  p <- add_scalebar(p, plot_df, length_um = 500, flip_y = flip_y, corner = "lr")
  
  out_file <- file.path(out_dir, paste0("segmentation_", make.names(im), ".png"))
  ggsave(out_file, plot = p, width = 10, height = 10, dpi = 800, limitsize = FALSE)
  
  rm(seg, plot_df, p); gc()
}
