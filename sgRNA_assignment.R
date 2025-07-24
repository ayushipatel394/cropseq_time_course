# sgRNA_assignment.R
# Author: Ayushi Patel
# Description: Assigns sgRNAs to cells from 10X scRNA-seq data and prepares a merged Seurat object

# Load required libraries
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(readr)
library(vegan)

# --------------------------------------------------------------------------------------------------
# Function to assign the most abundant sgRNA per cell
assign_gRNAs <- function(seurat_obj) {
  # Identify sgRNA features
  grna_list <- grep("sg", rownames(seurat_obj), value = TRUE)
  grna_list <- grep("-gene", grna_list, value = TRUE)

  if (length(grna_list) == 0) {
    stop("No gRNA features found matching 'sg' and '-gene' patterns.")
  }

  # Extract sgRNA expression matrix
  counts <- seurat_obj[["RNA"]]$counts
  grna_matrix <- counts[grna_list, ]

  if (nrow(grna_matrix) == 0) {
    stop("gRNA matrix is empty. Check gRNA feature names.")
  }

  # Assign most abundant sgRNA per cell
  transcriptome_assignments <- apply(grna_matrix, 2, function(x) {
    if (any(x > 0)) {
      max_val <- max(x)
      names(x)[which(x == max_val)][1]
    } else {
      NA
    }
  })

  # Create assignment dataframe
  transcriptome_df <- data.frame(
    gRNA = unname(transcriptome_assignments),
    cell_barcode = colnames(grna_matrix),
    stringsAsFactors = FALSE
  )

  # Clean and filter barcodes
  transcriptome_df <- transcriptome_df %>%
    separate(cell_barcode, "-", into = c("cellbarcode", "number"), remove = FALSE) %>%
    filter(as.numeric(number) < 12) %>%
    mutate(
      number = 1,
      cell_barcode = paste0(cellbarcode, "-", number),
      cell_barcode = gsub(" ", "", cell_barcode)
    )

  # Final assignment table
  assignment <- transcriptome_df %>%
    select(cell_barcode, gRNA) %>%
    separate(gRNA, "-", into = c("gene", "number", "random"), remove = FALSE)

  # Summary
  total_cells <- nrow(assignment)
  non_na_genes <- sum(!is.na(assignment$gene))
  message("Number of non-NA gene assignments: ", non_na_genes)
  message("Proportion of assigned cells: ", round(non_na_genes / total_cells, 3))

  return(assignment)
}

# --------------------------------------------------------------------------------------------------
# Import 10X matrices and create Seurat objects
data_c1t0 <- Read10X("/gpfs/home/ap6924/yanailab/projects/ap6924/ast_cropseq_poolc/KCL_poolc1_t0/outs/filtered_feature_bc_matrix")
data_c2t0 <- Read10X("/gpfs/home/ap6924/yanailab/projects/ap6924/ast_cropseq_poolc/KCL_poolc2_t0/outs/filtered_feature_bc_matrix")
data_c3t0 <- Read10X("/gpfs/home/ap6924/yanailab/projects/ap6924/ast_cropseq_poolc/KCL_poolc3_t0/outs/filtered_feature_bc_matrix")

c1t0 <- CreateSeuratObject(data_c1t0, project = "c1t0", assay = "RNA")
c2t0 <- CreateSeuratObject(data_c2t0, project = "c2t0", assay = "RNA")
c3t0 <- CreateSeuratObject(data_c3t0, project = "c3t0", assay = "RNA")

# --------------------------------------------------------------------------------------------------
# Assign sgRNAs and add to metadata
assignment_c1t0 <- assign_gRNAs(c1t0)
c1t0 <- AddMetaData(c1t0, metadata = assignment_c1t0)

assignment_c2t0 <- assign_gRNAs(c2t0)
c2t0 <- AddMetaData(c2t0, metadata = assignment_c2t0)

assignment_c3t0 <- assign_gRNAs(c3t0)
c3t0 <- AddMetaData(c3t0, metadata = assignment_c3t0)

# --------------------------------------------------------------------------------------------------
# Merge Seurat objects and clean up expression matrix
ct0 <- merge(x = c1t0, y = list(c2t0, c3t0))
ct0 <- JoinLayers(ct0)

# Remove sgRNA and Cas9 genes
counts <- LayerData(ct0, assay = "RNA", layer = "counts")
counts <- counts[-c(25240:26707), ]

# Keep only protein-coding genes
t3g <- read.table("/gpfs/home/ap6924/yanailab/projects/ap6924/transcriptional_modules/t2g.txt", header = TRUE)
t3g_proteincoding <- subset(t3g, transcript_biotype == "protein_coding")
counts <- counts[rownames(counts) %in% t3g_proteincoding$external_gene_name, ]

# Remove genes with 0 expression
counts2 <- counts[rowSums(counts) != 0, ]

# Subset Seurat object to filtered genes
ct0 <- subset(ct0, features = rownames(counts2))

# Save processed Seurat object
saveRDS(ct0, "poolct0_combine.rds")
