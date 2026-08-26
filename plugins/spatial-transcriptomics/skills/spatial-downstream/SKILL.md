---
name: spatial-downstream
description: Shared downstream analysis for all spatial platforms (Visium, Visium HD, Xenium, Atera) — neighborhood enrichment, CellChat cell-cell communication, and RNA velocity. Use after any platform branch has produced an h5ad with cell-type labels (classic Visium: deconvolution proportions; Visium HD / Xenium / Atera: cell-type annotations). Produces communication and trajectory results, then stops for review.
license: MIT
---

# Shared Downstream — Neighborhood, Communication, Velocity

## Goal

Run analyses that are **common across all platforms** after the platform branch produces an h5ad with cell-type labels. Granularity matches the input: spot-level for classic Visium, cell-level for Visium HD / Xenium / Atera.

## Prerequisites

- h5ad from any platform branch with:
  - `obsm['spatial']` (coordinates)
  - cell-type labels in `obs`:
    - Classic Visium: deconvolution proportions (dominant cell type per spot) from S4
    - Visium HD / Xenium / Atera: cell-type annotations (e.g., CellTypist / CellAssign)

## Steps

1. **Neighborhood enrichment**
   - `spatial_neighborhood_enrichment` (squidpy `nhood_enrichment`).
   - Granularity: spot-level (classic Visium) or cell-level (HD / Xenium / Atera) — implement per input granularity.
   - Plots: enrichment heatmap.
   - Review focus: are enriched co-localizations biologically plausible?

2. **CellChat cell-cell communication**
   - `infer_spatial_cell_communication` (CellChat).
   - Input labels: deconvolution proportions (classic Visium) or cell-type annotations (cell-level platforms).
   - Plots: CellChat network / heatmap plots.
   - Review focus: are inferred ligand-receptor interactions biologically plausible?

3. **RNA velocity (optional)**
   - Requires spliced/unspliced counts (e.g., from velocity-compatible pipelines).
   - **Check availability first**: classic Visium HD / Xenium / Atera outputs may NOT include spliced/unspliced layers. If absent, skip with explanation.
   - Plots: velocity stream on spatial coordinates / UMAP.

## Outputs

- `results/08_neighborhood/enrichment_zscore.csv` + heatmap
- `results/09_communication/CellChat_object.rds` + network plots
- `results/10_velocity/velocity_stream.png` (if applicable)

## Biological Interpretation

- Neighborhood: report top enriched co-localizations; relate to tissue architecture.
- CellChat: report top ligand-receptor pairs; relate to known biology (e.g., immune-tumor interactions).
- Velocity: report dominant trajectories if available.

## Stop for Review

Present interpretation using the template from the parent `spatial-transcriptomics` skill. Wait for `通过` / `调整` / `跳过` before final summary.

## Notes

- Granularity-aware: the same analysis runs at spot level for classic Visium and cell level for the other platforms — do NOT mix granularities.
- Velocity: only run if spliced/unspliced data exists; otherwise state clearly that it is skipped.
