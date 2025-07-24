# module_entropy_analysis_visualization.R
# --------------------------
# Compute Shannon entropy for modules across perturbations
# Author: [Your Name]
# Date: [Date]
# --------------------------

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggridges)
library(viridis)

# --------------------------
# 1. Setup and parameters
# --------------------------
meta <- combine@meta.data
meta$cell <- rownames(meta)
module_cols <- colnames(meta)[39:101]  # Adjust if needed

set.seed(123)
n_iter <- 100
n_cells_per_perturb <- 50

# --------------------------
# 2. Helper functions
# --------------------------

# Identify high-scoring cells in a module (mean + 1 SD)
get_high_scoring <- function(vec) {
  vec > (mean(vec, na.rm = TRUE) + sd(vec, na.rm = TRUE))
}

# Compute Shannon entropy (base-2)
compute_entropy <- function(labels) {
  probs <- table(labels) / length(labels)
  -sum(probs * log2(probs[probs > 0]))
}

# --------------------------
# 3. Main entropy computation
# --------------------------
all_entropy <- list()

for (i in seq_len(n_iter)) {
  sampled_cells <- meta %>%
    group_by(perturb) %>%
    sample_n(n_cells_per_perturb, replace = TRUE) %>%
    pull(cell)

  sampled_meta <- meta[sampled_cells, ]

  for (mod in module_cols) {
    high_cells <- get_high_scoring(sampled_meta[[mod]])
    if (sum(high_cells) < 10) next

    high_perturbs <- sampled_meta$perturb[high_cells]
    ent <- compute_entropy(high_perturbs)

    all_entropy[[length(all_entropy) + 1]] <- data.frame(
      Module = mod,
      Iteration = i,
      Entropy = ent
    )
  }
}

entropy_df <- bind_rows(all_entropy)
write.csv(entropy_df, "module_entropy_results.csv", row.names = FALSE)

# --------------------------
# 4. Entropy grouping and visualization
# --------------------------
# Group modules by mean entropy (quartiles)
mean_entropy <- entropy_df %>%
  group_by(Module) %>%
  summarise(MeanEntropy = mean(Entropy, na.rm = TRUE), .groups = "drop") %>%
  mutate(EntropyGroup = ntile(MeanEntropy, 4))

# Join group back to entropy dataframe
entropy_df <- left_join(entropy_df, mean_entropy, by = "Module")

# Density plots by entropy group
ggplot(entropy_df, aes(x = Entropy, fill = Module)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~EntropyGroup, ncol = 2, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Entropy Distribution by Module Group",
    x = "Shannon Entropy (bits)",
    y = "Density"
  ) +
  theme(legend.position = "none")

# --------------------------
# 5. Heatmap of mean entropy
# --------------------------
entropy_matrix <- mean_entropy %>%
  column_to_rownames("Module") %>%
  as.matrix() %>%
  t()  # Transpose for pheatmap

pheatmap(
  entropy_matrix,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  main = "Mean Entropy per Module"
)

# --------------------------
# 6. Boxplot of entropy per module
# --------------------------
# Summary stats for ordering
module_summary <- entropy_df %>%
  group_by(Module) %>%
  summarise(
    median_entropy = median(Entropy, na.rm = TRUE),
    mean_entropy = mean(Entropy, na.rm = TRUE),
    sd_entropy = sd(Entropy, na.rm = TRUE),
    q25 = quantile(Entropy, 0.25, na.rm = TRUE),
    q75 = quantile(Entropy, 0.75, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(median_entropy)

# Order modules by median
entropy_df$Module <- factor(entropy_df$Module, levels = module_summary$Module)

# Create boxplot
boxplot <- ggplot(entropy_df, aes(x = Module, y = Entropy, fill = Module)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(
    title = "Perturbation Specificity of Modules",
    subtitle = "Lower entropy = more perturbation-specific",
    x = "Module",
    y = "Shannon Entropy (bits)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

ggsave("module_entropy_boxplot.pdf", boxplot, width = 12, height = 7)

# --------------------------
# 7. Ridge plot for top 20 most specific modules
# --------------------------
top20_modules <- head(module_summary$Module, 20)
ridges_data <- entropy_df %>% filter(Module %in% top20_modules)

ridges_plot <- ggplot(ridges_data, aes(x = Entropy, y = Module, fill = Module)) +
  ggridges::geom_density_ridges(alpha = 0.7, scale = 0.9) +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(
    title = "Entropy Distribution of Top 20 Most Specific Modules",
    subtitle = "Lower entropy = more perturbation-specific",
    x = "Shannon Entropy (bits)",
    y = "Module"
  ) +
  theme(legend.position = "none")

ggsave("module_entropy_ridges.pdf", ridges_plot, width = 10, height = 8)
