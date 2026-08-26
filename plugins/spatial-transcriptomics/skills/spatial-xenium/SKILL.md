---
name: spatial-xenium
description: Xenium platform branch of the spatial transcriptomics workflow — load and validate the platform's cell-level matrix for downstream analysis. Use when the user's data is 10x Xenium output (cell_feature_matrix.h5, cell-level h5ad, or transcripts.csv + segmentation). Produces a cell-level h5ad, then stops for review.
license: MIT
---

# Xenium Branch — Cell-Level Loading & Validation

## Goal

Load Xenium data as a **cell-level** matrix. Xenium is an imaging-based in situ platform that outputs cell barcodes with platform-provided segmentation — **no re-segmentation and no deconvolution needed**.

## Prerequisites

- Xenium output, one of:
  - `cell_feature_matrix.h5` (standard Xenium output; cells × genes)
  - cell-level h5ad with `obsm['spatial']`
  - `transcripts.csv` (transcript-level) + segmentation boundaries
- Python: `scanpy`

## Steps

1. **Load cell matrix**
   - `sc.read_10x_h5(cell_feature_matrix.h5)` — cells × genes, sparse.
   - Attach coordinates: Xenium outputs have `pxl_col_in_fullres` / `pxl_row_in_fullres` or `obsm['spatial']`.

2. **Validate platform segmentation**
   - Confirm cell boundaries exist (nucleus or cell segmentation from Xenium pipeline).
   - Do NOT re-segment unless the user explicitly asks.

3. **QC & sanity check**
   - Report: number of cells, median genes/cell, median counts/cell.
   - Flag: very low counts (possible empty cells), very high counts (possible doublets), abnormal spatial distribution.

## Outputs

- `results/01_loading/<sample>_cells.h5ad` — cell-level AnnData
- `results/01_loading/QC_plots.png` — QC violin plots

## Biological Interpretation

- Report total cells and QC stats.
- Note tissue type and expected cell composition.
- Flag any platform QC warnings.

## Stop for Review

Present interpretation using the template from the parent `spatial-transcriptomics` skill. Wait for `通过` / `调整` / `跳过` before proceeding to shared downstream.

## Notes

- Xenium segmentation comes with the platform — treat as ground truth unless the user asks otherwise.
- Deconvolution is NOT needed for Xenium — cells are already resolved.
