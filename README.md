# scRNA-seq and scATAC-seq Integration Pipeline

A Seurat v5 pipeline that integrates separately measured scRNA-seq and scATAC-seq data, transfers cell type labels across modalities, and produces **paired multiome-style outputs** for downstream GRN inference.

## Why this exists

Most GRN inference tools (FigR, TRIPOD, Pando, etc.) expect paired multi-omic input where RNA and ATAC are measured in the same cells. This pipeline takes **unpaired** scRNA-seq and scATAC-seq data and produces 1:1 matched cell pairs via anchor-based integration — enabling these tools to run on separately collected datasets.

Built for benchmarking GRN methods across matched, mismatched, and unpaired scenarios.

## Input

Only **two CSV files** required (no fragments file needed):

| File | Format | Example row name |
|------|--------|-----------------|
| `rna_counts.csv` | Gene × cell matrix | `TP53` |
| `atac_counts.csv` | Peak × cell matrix | `chr1:115460-115953` or `chr1-115460-115953` |

## Output

| File | Description |
|------|-------------|
| `matched_rna_counts.csv` | Gene × N paired cells (raw counts) |
| `matched_atac_counts.csv` | Peak × N paired cells (same columns, same order) |
| `matched_cell_metadata.csv` | Barcode mapping, cell types, confidence scores |
| `rna_processed.rds` | Full RNA Seurat object |
| `atac_annotated.rds` | Full ATAC Seurat object with predicted labels |
| `coembedded.rds` | Merged object for joint visualization |

The two matched CSV files share identical column names (`Cell_0001`, `Cell_0002`, ...) so column N in RNA corresponds exactly to column N in ATAC.

## Quick start

```r
# Set your inputs and genome at the top of the script
RNA_CSV_PATH  <- "path/to/rna_counts.csv"
ATAC_CSV_PATH <- "path/to/atac_counts.csv"
GENOME        <- "hg38"    # or "mm10" for mouse

source("seurat5_atacseq_integration_from_csv.R")
```

## Configuration

All settings are at the top of the script:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `GENOME` | `"hg38"` | `"hg38"` for human, `"mm10"` for mouse (auto-configures EnsDb, mito pattern) |
| `CLUSTER_RESOLUTION` | `0.8` | Clustering resolution |
| `MIN_PREDICTION_SCORE` | `0.2` | Min confidence to keep a matched pair |

Peak separator is auto-detected (`chr1:start-end` vs `chr1-start-end`).

## How it works

1. **Load** CSV matrices into Seurat/ChromatinAssay objects
2. **Process RNA** — normalize, PCA, cluster, annotate cell types
3. **Process ATAC** — TF-IDF, LSI, cluster, attach gene annotations
4. **Bridge modalities** — overlap peaks with gene body + 2kb promoter to build a gene activity matrix (replaces `GeneActivity()`, no fragments file needed)
5. **Find anchors** — CCA on RNA expression vs ATAC gene activity
6. **Transfer labels** — predict RNA cell types onto ATAC cells
7. **Co-embed** — impute RNA into ATAC, merge, joint PCA/UMAP
8. **1:1 matching** — greedy best-score-first anchor pairing, each cell used once
9. **Export** — paired CSVs with raw counts ready for GRN tools

## Cell type annotation

The script uses cluster IDs as placeholder labels. Replace with your own annotations in Section 3:

```r
# Option A: Manual
markers <- FindAllMarkers(pbmc.rna, only.pos = TRUE)
new.ids <- c("CD14 Mono", "CD4 Naive", "B", ...)
names(new.ids) <- levels(pbmc.rna)
pbmc.rna <- RenameIdents(pbmc.rna, new.ids)

# Option B: Automated
pbmc.rna <- RunAzimuth(pbmc.rna, reference = "pbmcref")
pbmc.rna$cell_type <- pbmc.rna$predicted.celltype.l2
```

## Prerequisites

```r
install.packages(c("Seurat", "Signac", "ggplot2", "cowplot"))
BiocManager::install(c("EnsDb.Hsapiens.v86",    # human
                        "EnsDb.Mmusculus.v79",    # mouse
                        "GenomicRanges", "GenomeInfoDb"))
```

## Downstream compatibility

The paired outputs are compatible with GRN inference tools including:

- **Paired-cell tools**: FigR, TRIPOD, Pando, DIRECTNET
- **Unpaired-compatible tools**: SCENIC+, CellOracle, LINGER

## Limitations

- The 1:1 anchor matching is greedy, not globally optimal — consider this when interpreting GRN results
- Matching discards unmatched cells, reducing dataset size
- Gene activity from peak overlaps is an approximation of `GeneActivity()` with fragments
- Cell type annotation quality directly affects label transfer accuracy

## License

GNU
