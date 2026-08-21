---
name: spatial-s2-normalize-cluster
description: Stage 2 of the spatial transcriptomics workflow — normalize 10x Visium data and cluster spatial spots. Use when the user asks to normalize, log-transform, find highly variable genes, or cluster Visium spots (PCA/UMAP/Leiden). Produces UMAP and spatial cluster plots, then stops for review.
license: MIT
---

# S2 — Normalize + Cluster

## Goal

Normalize and log-transform the filtered data, then cluster spots with the standard PCA + UMAP + Leiden pipeline.

## Steps

1. **Normalize**: `normalize_visium(adata_path=..., output_path=..., plot_path=..., target_sum=..., n_top_genes=...)`
   - Defaults: `target_sum=None` (scanpy default), `n_top_genes=2000`.
   - Output: normalized h5ad + HVG plot.
2. **Cluster**: `cluster_spatial_data(adata_path=..., output_path=..., plot_path=..., n_pcs=..., resolution=...)`
   - Defaults: `n_pcs=30`, `resolution=1.0` (adjust if clusters look over/under-split).
   - Output: clustered h5ad + UMAP plot + spatial scatter colored by cluster.

## Outputs

- `results/03_normalization/<sample>_normalized.h5ad`
- `results/04_clustering/<sample>_clustered.h5ad`
- `results/04_clustering/<sample>_umap.png`
- `results/04_clustering/<sample>_spatial_clusters.png`

## Biological Interpretation

- Report the number of clusters found and their spatial distribution.
- Look at the spatial scatter: do clusters form contiguous spatial regions (tissue architecture) or scattered noise?
- If clusters are too fragmented / too coarse, suggest resolution adjustment.

## Stop for Review

Present interpretation using the template from the parent `spatial-transcriptomics` skill. Wait for `通过` / `调整` / `跳过` before proceeding to S3.
