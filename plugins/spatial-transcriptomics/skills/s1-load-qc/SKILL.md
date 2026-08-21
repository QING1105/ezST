---
name: spatial-s1-load-qc
description: Stage 1 of the spatial transcriptomics workflow — load 10x Visium data and QC-filter low-quality spots. Use when the user asks to load Visium data, initialize a spatial project, or filter spots by quality (counts, genes, mitochondrial fraction). Produces QC plots and a filtered AnnData, then stops for review.
license: MIT
---

# S1 — Load + QC

## Goal

Initialize the project, load Visium data into a standardized AnnData, and filter low-quality spots.

## Steps

1. **Init project**: `init_spatial_project(project_dir=...)` — creates the standardized directory layout.
2. **Load data**: `load_visium_data(counts_file=..., coordinates_file=..., sample_id=..., output_path=...)`
   - h5ad: pass the h5ad path as `counts_file` (coordinates already embedded).
   - 10x matrix: pass matrix/barcodes/features paths.
   - wide CSV: pass counts CSV + spot coordinates CSV.
   - Output: h5ad with `X` = raw counts, `obsm['spatial']`, `obs['barcode']`.
3. **Filter spots**: `filter_visium_spots(adata_path=..., output_path=..., plot_path=..., min_counts=..., min_genes=..., pct_mt=...)`
   - Defaults: `min_counts=500`, `min_genes=250`, `pct_mt=20` (adjust if user specifies).
   - Output: filtered h5ad + QC violin/spatial plots.

## Outputs

- `results/01_loading/<sample>_loaded.h5ad`
- `results/02_qc/<sample>_filtered.h5ad`
- `results/02_qc/<sample>_qc_*.png` (pre/post filter plots)

## Biological Interpretation

- Report spots retained: `filtered n / raw n` and the fraction removed.
- Look at spatial plots: is the removed area consistent with tissue edges / low-quality regions?
- State whether thresholds were reasonable or should be tightened/loosened.

## Stop for Review

Present interpretation using the template from the parent `spatial-transcriptomics` skill. Wait for `通过` / `调整` / `跳过` before proceeding to S2.
