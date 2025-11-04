# neighborhood_analysis.R
# Author: Erin McCaffrey 
# Last updated: 250813
# Overview: This script takes the clustered neighborhood data and performs several
# analysis. First it determines the neighborhood frequency per sample and visualizes
# as a stacked bar. Next it summarizies the neighborhood frequency per time
# point with statistical evaluation.

library(ggpubr)
library(dplyr)
library(ggplot2)
library(rstatix)

##...Load in data..##

setwd("/Volumes/T7 Shield/dang_lab_crypto_grans/Xenium/")
neighborhood_data <- read.csv('dang_xenium_neighborhoods_r50_k10.csv')
neighborhood_color_key <- read.csv('r50_k10_color-key.csv')

##..Get the prevalence of the clusters per tissue..##

neighborhood_summary <- neighborhood_data %>% 
  group_by(Samples) %>%
  count(as.factor(neighborhood_cluster), .drop=F) %>%
  mutate(freq_of_total = prop.table(n)) 

colnames(neighborhood_summary) <- c('sample','neighborhood_cluster','total',
                                    'frequency')

# append the group annotation and rename cluster column
neighborhood_summary <- neighborhood_summary %>%
  mutate(group = sub("_.*", "", sample))

##..Create a stacked bar across the samples..##

# reorder by descending median 
cluster_order <- neighborhood_summary %>%
  group_by(neighborhood_cluster) %>%
  summarise(median_freq = median(frequency), .groups = "drop") %>%
  arrange(desc(median_freq)) %>%
  pull(neighborhood_cluster)

neighborhood_summary <- neighborhood_summary %>%
  mutate(neighborhood_cluster = factor(neighborhood_cluster, levels = cluster_order))

# get the cluster colors
neighborhood_clusters <- levels(neighborhood_summary$neighborhood_cluster)
plot_colors <- droplevels(neighborhood_color_key[neighborhood_color_key$neighborhood_cluster %in% neighborhood_clusters,])
plot_colors$neighborhood_cluster <- factor(plot_colors$neighborhood_cluster, levels = neighborhood_clusters)
plot_colors <- plot_colors[order(plot_colors$neighborhood_cluster),]
color <- as.vector(plot_colors$hex)

# my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
#            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

ggplot(neighborhood_summary, aes(x = sample, y = frequency, fill = neighborhood_cluster)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = color) +
  labs(x = "sample", y = "neighborhood frequency", fill = "neighborhood cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##..Compare frequency across timepoints..##

# plot without p-values
ggplot(neighborhood_summary, aes(x = group, y = frequency, fill = neighborhood_cluster)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.6) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8) +  
  facet_wrap(~neighborhood_cluster, scales = "free_y", ncol = 5) +
  stat_summary(fun = median, geom = "line", aes(group = neighborhood_cluster),color = "black") +
  scale_fill_manual(values = color) +
  theme_bw() +
  labs(x = "timepoint", y = "neighborhood frequency",color = "neighborhood cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold")) +
  theme(legend.position = "None")

# add p-values
stat_tbl <- neighborhood_summary %>%
  group_by(neighborhood_cluster) %>%
  pairwise_wilcox_test(frequency ~ group, p.adjust.method = "none") %>%
  add_xy_position(x = "group") %>%
  mutate(y.position = y.position * 0.75)

ggplot(neighborhood_summary, aes(x = group, y = frequency, fill = neighborhood_cluster)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.6) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8) +  
  facet_wrap(~neighborhood_cluster, scales = "free_y", ncol = 5) +
  stat_summary(fun = median, geom = "line", aes(group = neighborhood_cluster),color = "black") +
  scale_fill_manual(values = color) +
  theme_bw() +
  labs(x = "timepoint", y = "neighborhood frequency",color = "neighborhood cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold")) +
  theme(legend.position = "None") + 
  stat_pvalue_manual(stat_tbl, label = "p.adj.signif", tip.length = 0, hide.ns = TRUE, step.increase = 0.005)


