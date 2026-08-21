---
name: spatial-s5-downstream
description: Stage 5 of the spatial transcriptomics workflow — neighborhood enrichment and cell-cell communication analysis. Use when the user asks to compute cell-type neighborhood enrichment (co-localization) or infer cell-cell communication with CellChat on Visium data. Produces enrichment heatmaps and CellChat network plots, then stops for review.
license: MIT
---

# S5 — Downstream (Neighborhood + Communication)

## Goal

Quantify cell-type co-localization with neighborhood enrichment, and infer spatially proximal cell-cell communication with CellChat v2.

## Prerequisites

- Deconvolved Visium h5ad from S4 with cell-type proportions (for both analyses).
- R + CellChat v2 for communication (see `install_r_packages_spatial.R`).

## Steps

1. **Neighborhood enrichment**: `spatial_neighborhood_enrichment(adata_path=..., cell_type_key=..., output_path=..., plot_path=..., n_neighs=...)`
   - Defaults: `n_neighs=6`.
   - Output: enrichment z-score table + heatmap.
2. **CellChat**: `infer_spatial_cell_communication(st_h5ad=..., cell_type_key=..., output_dir=..., plot_path=..., r_script_path=..., species=...)`
   - `species`: `"human"` (default) or `"mouse"` — selects the CellChatDB.
   - Output: CellChat object, communication summary, network/heatmap plots.

## Outputs

- `results/08_neighborhood/<sample>_enrichment.csv` + heatmap PNG
- `results/09_communication/<sample>_cellchat_*.png` (network, heatmap, bubble)
- `results/09_communication/<sample>_communication_summary.csv`

## Biological Interpretation

- Neighborhood: which cell-type pairs co-localize significantly (|z| > 2)? Does this match known tissue niches (e.g., tumor–macrophage, T cell–APC)?
- Communication: which ligand–receptor pairs dominate? Are they consistent with the tissue's known biology (e.g., immune checkpoint pairs in tumor)?

## Stop for Review

Present interpretation using the template from the parent `spatial-transcriptomics` skill. This is the final stage — after review, produce the overall workflow summary.
