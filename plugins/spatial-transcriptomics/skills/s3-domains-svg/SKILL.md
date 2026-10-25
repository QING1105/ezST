---
name: spatial-s3-domains-svg
description: Stage 3 of the spatial transcriptomics workflow — identify spatial domains and detect spatially variable genes. Use when the user asks to find spatial domains/niches, run Moran's I, or detect spatially variable genes (SVGs) in Visium data. Produces spatial domain maps and SVG rankings, then stops for review.
license: MIT
---

## 0.2.1 运行约定

表达聚类使用 `identify_spatial_domains(..., spatial_weight=0)`；加入空间邻接时显式设置 0 到 1 之间的权重，并报告 `uns['domain_clustering']`。CellChat 默认 `mode="expression"`；空间模式需要 `mode="spatial"`、`microns_per_coordinate`、`spatial_tolerance`，按样本尺度设置作用范围和距离缩放。输入优先使用保存的原始 counts。CellChat 输出实际为 `cellchat_communication.csv`、`cellchat_pathways.csv`、`cellchat.rds`、`cellchat_network.png` 与 `run_settings.json`。

StarDist 几何使用数值 `.npz`，旧 pickle 读取需显式 `allow_legacy_pickle=True`。HD 裁图通过 `--image`、`--cells`、`--geometry`、`--roi` 指定输入；全部坐标使用全分辨率图像像素。Proseg 基因名来自其导出的 `proseg_genes.csv.gz`，导入后核对 `var_names`。Python 与插件内运行文件通过 `python scripts/sync_plugin_runtime.py` 同步。



# S3 — Spatial Domains + SVG

## Goal

Cluster expression profiles, optionally blending a spatial neighborhood graph, and detect spatially variable genes with Moran's I. Interpret anatomical regions using marker and spatial evidence.

## Steps

1. **Spatial domains**: `identify_spatial_domains(adata_path=..., output_path=..., plot_path=..., resolution=..., n_neighs=..., spatial_weight=...)`
   - Defaults: `resolution=1.0`, `spatial_weight=0`; `n_neighs` is inferred from sample size. Report expression-only clustering when weight is zero.
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
