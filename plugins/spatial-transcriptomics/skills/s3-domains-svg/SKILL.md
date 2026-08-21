---
name: spatial-s3-domains-svg
description: Stage 3 of the spatial transcriptomics workflow — identify spatial domains and detect spatially variable genes. Use when the user asks to find spatial domains/niches, run Moran's I, or detect spatially variable genes (SVGs) in Visium data. Produces spatial domain maps and SVG rankings, then stops for review.
license: MIT
---

# S3 — Spatial Domains + SVG

## Goal

Identify spatial domains (tissue architecture regions) using a spatial neighborhood graph, and detect spatially variable genes with Moran's I.

## Steps

1. **Spatial domains**: `identify_spatial_domains(adata_path=..., output_path=..., plot_path=..., resolution=..., n_neighs=...)`
   - Defaults: `resolution=1.0`, `n_neighs=6` (spatial neighbor count).
   - Output: domain assignments + spatial domain map.
2. **SVG detection**: `find_spatially_variable_genes(adata_path=..., output_path=..., plot_path=..., n_top_genes=..., n_jobs=..., genes_to_plot=...)`
   - Defaults: `n_top_genes=100`, `n_jobs=1`.
   - `genes_to_plot`: optional list of specific genes to visualize on spatial coordinates (useful for marker validation, e.g., `EPCAM`, `CD3D`, `VIM`).
   - Output: SVG ranked CSV + spatial plots of top genes.

## Outputs

- `results/05_domains/<sample>_domains.h5ad`
- `results/05_domains/<sample>_domains.png`
- `results/06_svg/<sample>_svg_results.csv`
- `results/06_svg/<sample>_svg_*.png`

## Biological Interpretation

- Report the number of spatial domains and their spatial layout.
- Interpret domains against expected tissue architecture (e.g., tumor core, stroma, immune infiltrate, normal epithelium).
- Check top SVGs: are they known markers for the tissue type? Do their spatial patterns align with domains?

## Stop for Review

Present interpretation using the template from the parent `spatial-transcriptomics` skill. Wait for `通过` / `调整` / `跳过` before proceeding to S4.
