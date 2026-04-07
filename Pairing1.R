#!/usr/bin/env Rscript
# =============================================================================
# Integrating scRNA-seq and scATAC-seq data (Seurat v5)
# Works from two CSV files ONLY — no fragments file needed
# =============================================================================
#
# INPUT (only 2 files):
#   1. rna_counts.csv   — Gene-by-cell matrix (rows = genes, cols = barcodes)
#   2. atac_counts.csv  — Peak-by-cell matrix (rows = peaks, cols = barcodes)
#                          Peak names MUST be "chr-start-end" (dash-separated)
#
# NOTE: The original Seurat vignette uses GeneActivity() which requires a
#   fragments.tsv.gz file. This script REPLACES that step with a manual
#   peak-to-gene overlap approach using only the peak count matrix.
#   It overlaps peaks with gene promoter + body regions and sums peak counts
#   per gene to create a gene activity matrix.
#
# GENOME: Assumes human hg38. For mouse, swap:
#   - EnsDb.Hsapiens.v86 -> EnsDb.Mmusculus.v79
#   - genome "hg38" -> "mm10"
#
# Prerequisites:
#   install.packages(c("Seurat", "Signac", "ggplot2", "cowplot"))
#   BiocManager::install(c("EnsDb.Hsapiens.v86", "GenomicRanges", "GenomeInfoDb"))
# =============================================================================


library(Seurat)
library(Signac)
#library(EnsDb.Hsapiens.v86)
#library(EnsDb.Mmusculus.v79)
library(GenomicRanges)
library(GenomeInfoDb)
library(Matrix)
library(ggplot2)
library(cowplot)

# ========================== USER SETTINGS ====================================
GENOME <- "hg38"   # change to "mm10" for mouse

# ========================== AUTO-CONFIGURE BY GENOME =========================
if (GENOME == "mm10") {
  library(EnsDb.Mmusculus.v79)
  ENSDB       <- EnsDb.Mmusculus.v79
  MT_PATTERN  <- "^mt-"
} else {
  library(EnsDb.Hsapiens.v86)
  ENSDB       <- EnsDb.Hsapiens.v86
  MT_PATTERN  <- "^MT-"
}

# ========================== USER SETTINGS ====================================
RNA_CSV_PATH  <- "/Path/to/subsampled_RNA_Macrophase_buffer1_filtered.csv"
ATAC_CSV_PATH <- "/Path/to/subsampled_ATAC_K562_human_filtered.csv"


# QC thresholds
#RNA_MIN_FEATURES  <- 200
#RNA_MAX_FEATURES  <- 5000
#RNA_MAX_MT_PCT    <- 20
#ATAC_MIN_FEATURES <- 500
#ATAC_MAX_FEATURES <- 50000

CLUSTER_RESOLUTION <- 0.8

# Promoter extension upstream of TSS (bp) for gene activity calculation
PROMOTER_UPSTREAM <- 2000
PROMOTER_DOWNSTREAM <- 0

# Minimum prediction score to keep a matched pair
MIN_PREDICTION_SCORE <- 0.2
# =============================================================================


# =============================================================================
# SECTION 1: Load CSVs and create Seurat objects
# =============================================================================

cat("=== Loading RNA count matrix ===\n")
rna_counts <- read.csv(RNA_CSV_PATH, row.names = 1, check.names = FALSE)
rna_counts <- as(as.matrix(rna_counts), "dgCMatrix")
cat("  Dimensions:", nrow(rna_counts), "genes x", ncol(rna_counts), "cells\n")

pbmc.rna <- CreateSeuratObject(
  counts  = rna_counts,
  assay   = "RNA",
  project = "RNA"
)
rm(rna_counts); gc()

cat("\n=== Loading ATAC count matrix ===\n")
atac_counts <- read.csv(ATAC_CSV_PATH, row.names = 1, check.names = FALSE)
atac_counts <- as(as.matrix(atac_counts), "dgCMatrix")
cat("  Dimensions:", nrow(atac_counts), "peaks x", ncol(atac_counts), "cells\n")

first_peak <- rownames(atac_counts)[1]
if (grepl(":", first_peak)) {
  peak_sep <- c(":", "-")
} else {
  peak_sep <- c("-", "-")
}

atac_assay <- CreateChromatinAssay(
  counts       = atac_counts,
  sep          = peak_sep,
  genome       = GENOME,
  min.cells    = 0,
  min.features = 0
)

pbmc.atac <- CreateSeuratObject(
  counts  = atac_assay,
  assay   = "peaks",
  project = "ATAC"
)
# Keep a copy of the raw ATAC counts for later export
atac_counts_backup <- atac_counts
rm(atac_counts); gc()


# =============================================================================
# SECTION 2: Quality control
# =============================================================================

#cat("\n=== QC filtering ===\n")

##pbmc.rna[["percent.mt"]] <- PercentageFeatureSet(pbmc.rna, pattern = "^MT-")
#cat("  RNA cells before QC:", ncol(pbmc.rna), "\n")
#pbmc.rna <- subset(pbmc.rna,
#                    nFeature_RNA > RNA_MIN_FEATURES &
#                    nFeature_RNA < RNA_MAX_FEATURES &
#                    percent.mt   < RNA_MAX_MT_PCT)
#cat("  RNA cells after QC: ", ncol(pbmc.rna), "\n")
#
#cat("  ATAC cells before QC:", ncol(pbmc.atac), "\n")
#pbmc.atac <- subset(pbmc.atac,
#                     nFeature_peaks > ATAC_MIN_FEATURES &
#                     nFeature_peaks < ATAC_MAX_FEATURES)
#cat("  ATAC cells after QC: ", ncol(pbmc.atac), "\n")
#
# Update backup to match QC-filtered cells
atac_counts_backup <- atac_counts_backup[, colnames(pbmc.atac)]


# =============================================================================
# SECTION 3: Process RNA
# =============================================================================

cat("\n=== Processing RNA ===\n")
pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna, nfeatures = 2000)
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna, verbose = FALSE)
pbmc.rna <- FindNeighbors(pbmc.rna, dims = 1:30)
pbmc.rna <- FindClusters(pbmc.rna, resolution = CLUSTER_RESOLUTION)
pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)

# ---- CELL TYPE ANNOTATION ----
# REPLACE with your own annotations!
pbmc.rna$cell_type <- paste0("Cluster_", Idents(pbmc.rna))
cat("  Cell types:", length(unique(pbmc.rna$cell_type)), "groups\n")


# =============================================================================
# SECTION 4: Process ATAC
# =============================================================================

cat("\n=== Processing ATAC ===\n")

# Add gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = ENSDB)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- GENOME
Annotation(pbmc.atac) <- annotations

pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")
pbmc.atac <- RunSVD(pbmc.atac)
pbmc.atac <- FindNeighbors(pbmc.atac, reduction = "lsi", dims = 2:30)
pbmc.atac <- FindClusters(pbmc.atac, resolution = CLUSTER_RESOLUTION)
pbmc.atac <- RunUMAP(pbmc.atac,
                      reduction      = "lsi",
                      dims           = 2:30,
                      reduction.name = "umap.atac",
                      reduction.key  = "atacUMAP_")


# =============================================================================
# SECTION 5: Build gene activity matrix from peak overlaps
#             (REPLACES GeneActivity() — no fragments file needed)
# =============================================================================

cat("\n=== Building gene activity matrix from peak-gene overlaps ===\n")

# Step 1: Get gene coordinates and extend to include promoter
gene_coords <- genes(ENSDB, columns = "gene_name")
seqlevelsStyle(gene_coords) <- "UCSC"
gene_coords <- keepStandardChromosomes(gene_coords, pruning.mode = "coarse")

# Extend gene body upstream to include promoter region
gene_coords_extended <- gene_coords
start(gene_coords_extended) <- ifelse(
  strand(gene_coords_extended) == "+",
  start(gene_coords_extended) - PROMOTER_UPSTREAM,
  start(gene_coords_extended) - PROMOTER_DOWNSTREAM
)
end(gene_coords_extended) <- ifelse(
  strand(gene_coords_extended) == "+",
  end(gene_coords_extended) + PROMOTER_DOWNSTREAM,
  end(gene_coords_extended) + PROMOTER_UPSTREAM
)
# Ensure no negative starts
start(gene_coords_extended) <- pmax(start(gene_coords_extended), 1)

# Step 2: Get peak coordinates from the ATAC assay
peak_names <- rownames(pbmc.atac)
peak_split <- strsplit(peak_names, "-")

peak_gr <- GRanges(
  seqnames = sapply(peak_split, function(x) paste(x[1:(length(x)-2)], collapse = "-")),
  ranges   = IRanges(
    start = as.integer(sapply(peak_split, function(x) x[length(x)-1])),
    end   = as.integer(sapply(peak_split, function(x) x[length(x)]))
  )
)
names(peak_gr) <- peak_names

# Step 3: Find overlaps between peaks and gene regions
# Focus only on variable genes from RNA to match the original vignette
var_genes <- VariableFeatures(pbmc.rna)
gene_idx  <- which(gene_coords_extended$gene_name %in% var_genes)
gene_sub  <- gene_coords_extended[gene_idx]

# Remove duplicate gene names (keep first occurrence)
gene_sub <- gene_sub[!duplicated(gene_sub$gene_name)]

cat("  Variable genes with coordinates:", length(gene_sub), "\n")

hits <- findOverlaps(query = peak_gr, subject = gene_sub)

cat("  Peak-gene overlaps found:", length(hits), "\n")

# Step 4: Build gene activity matrix by summing peak counts per gene
peak_counts <- GetAssayData(pbmc.atac, assay = "peaks", layer = "counts")

# Create a sparse gene-by-cell matrix
gene_names_out <- gene_sub$gene_name
n_genes <- length(gene_names_out)
n_cells <- ncol(peak_counts)

# Build the matrix using overlap indices
# For each gene, sum the counts of all overlapping peaks
gene_activity <- sparseMatrix(
  i = integer(0), j = integer(0), x = numeric(0),
  dims = c(n_genes, n_cells),
  dimnames = list(gene_names_out, colnames(peak_counts))
)

cat("  Summing peak counts per gene...\n")

# Group overlaps by gene (subject hits)
subject_hits <- subjectHits(hits)
query_hits   <- queryHits(hits)

for (g in seq_len(n_genes)) {
  peak_indices <- query_hits[subject_hits == g]
  if (length(peak_indices) > 0) {
    if (length(peak_indices) == 1) {
      gene_activity[g, ] <- peak_counts[peak_indices, ]
    } else {
      gene_activity[g, ] <- colSums(peak_counts[peak_indices, , drop = FALSE])
    }
  }
}

# Remove genes with zero total activity
gene_sums <- rowSums(gene_activity)
gene_activity <- gene_activity[gene_sums > 0, ]
cat("  Gene activity matrix:", nrow(gene_activity), "genes x", ncol(gene_activity), "cells\n")

# Step 5: Add as a new assay to the ATAC object
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene_activity)
DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac, features = rownames(pbmc.atac))

cat("  Gene activity matrix built successfully (no fragments file needed)\n")


# =============================================================================
# SECTION 6: Find anchors and transfer labels
# =============================================================================

cat("\n=== Finding transfer anchors (CCA) ===\n")

# Use only genes present in both RNA variable features and the gene activity matrix
common_features <- intersect(VariableFeatures(pbmc.rna), rownames(gene_activity))
cat("  Common features for CCA:", length(common_features), "\n")

transfer.anchors <- FindTransferAnchors(
  reference       = pbmc.rna,
  query           = pbmc.atac,
  features        = common_features,
  reference.assay = "RNA",
  query.assay     = "ACTIVITY",
  reduction       = "cca"
)

cat("\n=== Transferring cell type labels ===\n")
celltype.predictions <- TransferData(
  anchorset        = transfer.anchors,
  refdata          = pbmc.rna$cell_type,
  weight.reduction = pbmc.atac[["lsi"]],
  dims             = 2:30
)

pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
pbmc.atac$cell_type <- pbmc.atac$predicted.id


# =============================================================================
# SECTION 7: Visualize
# =============================================================================

cat("\n=== Generating plots ===\n")

p1 <- DimPlot(pbmc.rna, group.by = "cell_type", label = TRUE) +
  NoLegend() + ggtitle("RNA (original labels)")
p2 <- DimPlot(pbmc.atac, group.by = "cell_type", label = TRUE) +
  NoLegend() + ggtitle("ATAC (predicted labels)")
print(p1 + p2)

p3 <- ggplot(pbmc.atac[[]], aes(x = prediction.score.max)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  theme_cowplot() +
  xlab("Max Prediction Score") + ylab("Number of cells") +
  ggtitle("Label transfer confidence")
print(p3)


# =============================================================================
# SECTION 8: Co-embedding
# =============================================================================

cat("\n=== Co-embedding RNA + ATAC ===\n")

genes.use <- VariableFeatures(pbmc.rna)
refdata   <- GetAssayData(pbmc.rna, assay = "RNA", layer = "data")[genes.use, ]

imputation <- TransferData(
  anchorset        = transfer.anchors,
  refdata          = refdata,
  weight.reduction = pbmc.atac[["lsi"]],
  dims             = 2:30
)

pbmc.atac[["RNA"]] <- imputation

coembed <- merge(x = pbmc.rna, y = pbmc.atac)
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

p4 <- DimPlot(coembed, group.by = "orig.ident") + ggtitle("By modality")
p5 <- DimPlot(coembed, group.by = "cell_type", label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("By cell type")
print(p4 + p5)


# =============================================================================
# SECTION 9: 1:1 cell matching from anchors
# =============================================================================

cat("\n=== Building 1:1 cell pairs from anchors ===\n")

anchor_df <- as.data.frame(transfer.anchors@anchors)
rna_cells  <- colnames(pbmc.rna)
atac_cells <- colnames(pbmc.atac)
anchor_df$rna_barcode  <- rna_cells[anchor_df$cell1]
anchor_df$atac_barcode <- atac_cells[anchor_df$cell2]

cat("  Total anchor pairs:", nrow(anchor_df), "\n")

# Greedy 1:1 matching (best score first)
anchor_df <- anchor_df[order(-anchor_df$score), ]

used_rna  <- character(0)
used_atac <- character(0)
matched_pairs <- data.frame(
  rna_barcode  = character(0),
  atac_barcode = character(0),
  anchor_score = numeric(0),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(anchor_df))) {
  rna_bc  <- anchor_df$rna_barcode[i]
  atac_bc <- anchor_df$atac_barcode[i]
  if (!(rna_bc %in% used_rna) && !(atac_bc %in% used_atac)) {
    matched_pairs <- rbind(matched_pairs, data.frame(
      rna_barcode  = rna_bc,
      atac_barcode = atac_bc,
      anchor_score = anchor_df$score[i],
      stringsAsFactors = FALSE
    ))
    used_rna  <- c(used_rna, rna_bc)
    used_atac <- c(used_atac, atac_bc)
  }
}

cat("  1:1 pairs before filtering:", nrow(matched_pairs), "\n")

# Add metadata
matched_pairs$prediction_score <- pbmc.atac$prediction.score.max[matched_pairs$atac_barcode]
matched_pairs$cell_type        <- pbmc.atac$cell_type[matched_pairs$atac_barcode]
matched_pairs$rna_cell_type    <- pbmc.rna$cell_type[matched_pairs$rna_barcode]
matched_pairs$labels_agree     <- matched_pairs$rna_cell_type == matched_pairs$cell_type

# Filter by prediction score
matched_pairs <- matched_pairs[matched_pairs$prediction_score >= MIN_PREDICTION_SCORE, ]
cat("  1:1 pairs after filtering (score >=", MIN_PREDICTION_SCORE, "):", nrow(matched_pairs), "\n")
cat("  Label agreement:", round(mean(matched_pairs$labels_agree) * 100, 1), "%\n")


# =============================================================================
# SECTION 10: Export matched multiome-style CSVs
# =============================================================================

cat("\n=== Exporting matched CSVs ===\n")

n_matched   <- nrow(matched_pairs)
unified_ids <- paste0("Cell_", sprintf("%04d", seq_len(n_matched)))

# RNA counts for matched cells
DefaultAssay(pbmc.atac) <- "peaks"

rna_matched  <- GetAssayData(pbmc.rna, assay = "RNA", layer = "counts")[, matched_pairs$rna_barcode]
atac_matched <- atac_counts_backup[, matched_pairs$atac_barcode]

colnames(rna_matched)  <- unified_ids
colnames(atac_matched) <- unified_ids


# Convert to dense and write CSV
cat("  Writing matched_rna_counts.csv  (", nrow(rna_matched), " genes x ", n_matched, " cells)...\n")
write.csv(as.matrix(rna_matched),  "/Path/to/matched_MacS1_RNA_K562_ATAC_rna_counts.csv")

cat("  Writing matched_atac_counts.csv (", nrow(atac_matched), " peaks x ", n_matched, " cells)...\n")
write.csv(as.matrix(atac_matched), "/Output/path/matched_MacS1_RNA_K562_ATAC_atac_counts.csv")



matched_meta <- data.frame(
  unified_cell_id  = unified_ids,
  rna_barcode      = matched_pairs$rna_barcode,
  atac_barcode     = matched_pairs$atac_barcode,
  cell_type        = matched_pairs$cell_type,
  rna_cell_type    = matched_pairs$rna_cell_type,
  labels_agree     = matched_pairs$labels_agree,
  anchor_score     = matched_pairs$anchor_score,
  prediction_score = matched_pairs$prediction_score,
  stringsAsFactors = FALSE
)
write.csv(matched_meta, "/Output/path/matched_MacS1_RNA_K562_ATAC_cell_metadata.csv", row.names = FALSE)


# =============================================================================
# SECTION 11: Save objects
# =============================================================================

#cat("\n=== Saving Seurat objects ===\n")
#saveRDS(pbmc.rna,  "rna_processed.rds")
#saveRDS(pbmc.atac, "atac_annotated.rds")
#saveRDS(coembed,   "coembedded.rds")


# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat("========================================================\n")
cat("  PIPELINE COMPLETE\n")
cat("========================================================\n")
cat("\n")
cat("  Input: 2 CSV files only (no fragments file needed)\n")
cat("    RNA cells:       ", ncol(pbmc.rna), "\n")
cat("    ATAC cells:      ", ncol(pbmc.atac), "\n")
cat("\n")
cat("  Matched output (same cells, same order):\n")
cat("    Paired cells:    ", n_matched, "\n")
cat("    matched_rna_counts.csv\n")
cat("    matched_atac_counts.csv\n")
cat("    matched_cell_metadata.csv\n")
cat("\n")
cat("  Seurat objects:\n")
cat("    rna_processed.rds\n")
cat("    atac_annotated.rds\n")
cat("    coembedded.rds\n")
cat("========================================================\n")
