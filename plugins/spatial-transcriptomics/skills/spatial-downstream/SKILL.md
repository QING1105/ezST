---
name: spatial-downstream
description: Shared downstream analysis for all spatial platforms (Visium, Visium HD, Xenium, Atera) — neighborhood enrichment, CellChat cell-cell communication, and RNA velocity. Use after any platform branch has produced an h5ad with cell-type labels (classic Visium: deconvolution proportions; Visium HD / Xenium / Atera: cell-type annotations). Produces communication and trajectory results, then stops for review.
license: MIT
---

## 0.2.1 运行约定

表达聚类使用 `identify_spatial_domains(..., spatial_weight=0)`；加入空间邻接时显式设置 0 到 1 之间的权重，并报告 `uns['domain_clustering']`。CellChat 默认 `mode="expression"`；空间模式需要 `mode="spatial"`、`microns_per_coordinate`、`spatial_tolerance`，按样本尺度设置作用范围和距离缩放。输入优先使用保存的原始 counts。CellChat 输出实际为 `cellchat_communication.csv`、`cellchat_pathways.csv`、`cellchat.rds`、`cellchat_network.png` 与 `run_settings.json`。

StarDist 几何使用数值 `.npz`，旧 pickle 读取需显式 `allow_legacy_pickle=True`。HD 裁图通过 `--image`、`--cells`、`--geometry`、`--roi` 指定输入；该 HD 裁图接口使用全分辨率图像像素坐标；Atera 数据保留其微米坐标。Proseg 基因名来自其导出的 `proseg_genes.csv.gz`，导入后核对 `var_names`。Python 与插件内运行文件通过 `python scripts/sync_plugin_runtime.py` 同步。



# Shared Downstream — Neighborhood, Communication, Velocity

## Goal

Run analyses that are **common across all platforms** after the platform branch produces an h5ad with cell-type labels. Granularity matches the input: spot-level for classic Visium, cell-level for Visium HD / Xenium / Atera.

## Prerequisites

- h5ad from any platform branch with:
  - `obsm['spatial']` (coordinates)
  - cell-type labels in `obs`:
    - Classic Visium: deconvolution proportions (dominant cell type per spot) from S4
    - Visium HD / Xenium / Atera: cell-type annotations (e.g., CellTypist / CellAssign)

## Atera input

For native Atera, use the reviewed cell-level h5ad from `spatial-atera`, including `layers["counts"]`, `obs["cell_type"]` and `obsm["spatial"]` in micrometres. Run separately for each ROI. Spatial CellChat uses `microns_per_coordinate=1`; report its actual distance parameters. Choose the annotation field and cell groups from the current sample, and record the selection used for downstream analysis.

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
- `cellchat.rds`, `cellchat_communication.csv`, `cellchat_pathways.csv`, `cellchat_network.png`, `run_settings.json` in the selected communication output directory
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
