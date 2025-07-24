# hotspot_module_analysis.R
# --------------------------
# Hypergeometric analysis between Hotspot modules and published transcriptional modules to annotate modules
# Author: Ayushi Patel
# --------------------------

# Load libraries
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# --------------------------
# 1. Load Hotspot module genes
# --------------------------
load_hotspot_modules <- function(directory) {
  csv_files <- list.files(path = directory, pattern = "module.*\\.csv$", full.names = TRUE)
  csv_data <- lapply(csv_files, read.csv)
  names(csv_data) <- gsub(".csv$", "", basename(csv_files))
  module_list <- lapply(csv_data, function(df) df$Gene)
  return(module_list)
}

hotspot_dir <- "/gpfs/data/yanailab/projects/ap6924/ast_cropseq_poolb/combine_control_analysis"
control_module_list <- load_hotspot_modules(hotspot_dir)
saveRDS(control_module_list, "../consensus_hotspot/control_module_list.RDS")

# --------------------------
# 2. Load published transcriptional modules
# --------------------------
db_mouse <- readRDS("../../transcriptional_modules/db_module_list_mouse.Rdata")
gavish_mouse <- readRDS("../../transcriptional_modules/gavish_modules_mouse.RDS")
tt_mouse <- readRDS("../../transcriptional_modules/ttmodules_mouse.RDS")

# --------------------------
# 3. Hypergeometric test function
# --------------------------
perform_hypergeometric_test <- function(list_genes, table_genes, total_genes = 3000) {
  q <- length(intersect(list_genes, table_genes))
  m <- length(table_genes)
  k <- length(list_genes)
  n <- total_genes - m
  phyper(q, m, n, k, lower.tail = FALSE)
}

# --------------------------
# 4. Run pairwise tests between module sets
# --------------------------
run_enrichment_tests <- function(reference_modules, test_modules, ref_label) {
  results <- expand.grid(List = names(reference_modules), Table = names(test_modules), stringsAsFactors = FALSE)
  results$P_value <- mapply(
    function(l, t) perform_hypergeometric_test(reference_modules[[l]], test_modules[[t]]),
    results$List, results$Table
  )
  results$Source <- ref_label
  return(results)
}

results_db <- run_enrichment_tests(db_mouse, control_module_list, "DB")
results_ag <- run_enrichment_tests(gavish_mouse, control_module_list, "AG")
results_tt <- run_enrichment_tests(tt_mouse, control_module_list, "TT")

# --------------------------
# 5. Combine results and plot heatmap
# --------------------------
all_results <- bind_rows(results_db, results_ag, results_tt)

# Convert to matrix form for heatmap
heatmap_matrix <- all_results %>%
  select(Source, List, Table, P_value) %>%
  unite(Module, Source, List, sep = "_") %>%
  pivot_wider(names_from = Module, values_from = P_value) %>%
  column_to_rownames("Table") %>%
  as.matrix()

# Remove unwanted rows (if applicable)
heatmap_matrix <- heatmap_matrix[!rownames(heatmap_matrix) %in% c("module6", "module7", "module17"), ]

# Plot heatmap
pdf("controlmodules_hypergeometric_p-value_rdbu.pdf", width = 18, height = 9)
pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(rev(brewer.pal(9, "Reds")))(50),
  breaks = seq(0, 0.05, length.out = 51),
  main = "p-value of hypergeometric test (0 to 0.05)",
  show_rownames = TRUE,
  show_colnames = TRUE
)
dev.off()

# --------------------------
# 6. Annotate modules with significant overlaps
# --------------------------
annotate_significant_modules <- function(results_df, alpha = 0.05) {
  significant <- results_df %>%
    filter(P_value < alpha) %>%
    group_by(Table) %>%
    summarise(Annotated = paste(List, collapse = "_"), .groups = "drop")

  # Add annotation to full list
  annotated_modules <- data.frame(Module = names(control_module_list)) %>%
    left_join(significant, by = c("Module" = "Table")) %>%
    mutate(Annotated = replace_na(Annotated, "Unannotated"))

  return(annotated_modules)
}

annotations_db <- annotate_significant_modules(results_db)
annotations_ag <- annotate_significant_modules(results_ag)
annotations_tt <- annotate_significant_modules(results_tt)

# Merge all annotations
all_annotations <- reduce(list(annotations_db, annotations_ag, annotations_tt), full_join, by = "Module") %>%
  unite(Final_Annotation, starts_with("Annotated"), sep = "_", na.rm = TRUE)

# Save annotations
write.csv(all_annotations, "control_module_annotations.csv", row.names = FALSE)
