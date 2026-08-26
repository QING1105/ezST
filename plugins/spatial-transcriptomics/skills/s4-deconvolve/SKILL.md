---
name: spatial-s4-deconvolve
description: Stage 4 of the spatial transcriptomics workflow — deconvolve Visium spots into cell-type proportions. Use when the user asks to deconvolve spatial spots, estimate cell-type composition, or run DestVI / SPOTlight with a scRNA-seq reference. Produces cell-type proportion maps, then stops for review.
license: MIT
---

# S4 — Deconvolution

## Goal

Estimate cell-type proportions for each Visium spot using a reference scRNA-seq dataset with cell-type labels.

## Prerequisites

- Reference scRNA-seq h5ad with cell-type labels in `obs` (e.g., column `level_3`).
- Clustered Visium h5ad from S2 (or raw loaded h5ad — DestVI works on raw counts).

## Steps

1. **DestVI** (preferred): `deconvolve_spatial_destvi(ref_h5ad=..., st_h5ad=..., cell_type_key=..., output_dir=..., plot_path=..., n_latent=..., max_epochs=..., destvi_max_epochs=...)`
   - Defaults: `n_latent=30`, `max_epochs=400` (CondSCVI), `destvi_max_epochs=500`.
   - Requires `scvi-tools` installed.
2. **SPOTlight** (alternative): `deconvolve_spatial_spotlight(ref_h5ad=..., st_h5ad=..., cell_type_key=..., output_dir=..., plot_path=..., r_script_path=...)`
   - Requires R + SPOTlight package (see `install_r_packages_spatial.R`).

## Outputs

- `results/07_deconvolution/<sample>_proportions.csv` (spots × cell types)
- `results/07_deconvolution/<sample>_proportions_*.png` (spatial proportion maps, stacked bar)

## Biological Interpretation

- Report dominant cell types and their spatial localization.
- Cross-check with tissue type: e.g., gastric cancer should show epithelial/tumor cells in tumor regions, immune cells at infiltrate margins.
- Flag any surprising dominant cell type (e.g., >50% of a cell type not expected in the tissue).

## Stop for Review

Present interpretation using the template from the parent `spatial-transcriptomics` skill. Wait for `通过` / `调整` / `跳过` before proceeding to S5.
