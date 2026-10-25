---
name: spatial-s5-downstream
description: Stage 5 of the spatial transcriptomics workflow — neighborhood enrichment and cell-cell communication analysis. Use when the user asks to compute cell-type neighborhood enrichment (co-localization) or infer cell-cell communication with CellChat on Visium data. Produces enrichment heatmaps and CellChat network plots, then stops for review.
license: MIT
---

## 0.2.1 运行约定

表达聚类使用 `identify_spatial_domains(..., spatial_weight=0)`；加入空间邻接时显式设置 0 到 1 之间的权重，并报告 `uns['domain_clustering']`。CellChat 默认 `mode="expression"`；空间模式需要 `mode="spatial"`、`microns_per_coordinate`、`spatial_tolerance`，按样本尺度设置作用范围和距离缩放。输入优先使用保存的原始 counts。CellChat 输出实际为 `cellchat_communication.csv`、`cellchat_pathways.csv`、`cellchat.rds`、`cellchat_network.png` 与 `run_settings.json`。

StarDist 几何使用数值 `.npz`，旧 pickle 读取需显式 `allow_legacy_pickle=True`。HD 裁图通过 `--image`、`--cells`、`--geometry`、`--roi` 指定输入；全部坐标使用全分辨率图像像素。Proseg 基因名来自其导出的 `proseg_genes.csv.gz`，导入后核对 `var_names`。Python 与插件内运行文件通过 `python scripts/sync_plugin_runtime.py` 同步。



# S5 — Downstream (Neighborhood + Communication)

## Goal

Quantify cell-type co-localization with neighborhood enrichment, and infer expression-based or explicitly distance-constrained cell-cell communication with CellChat v2.

## Prerequisites

- Deconvolved Visium h5ad from S4 with cell-type proportions (for both analyses).
- R + CellChat v2 for communication (see `install/install_r_packages.R`).

## Steps

1. **Neighborhood enrichment**: `spatial_neighborhood_enrichment(adata_path=..., cell_type_key=..., output_path=..., plot_path=..., n_neighs=...)`
   - Defaults: `n_neighs=6`.
   - Output: enrichment z-score table + heatmap.
2. **CellChat**: `infer_spatial_cell_communication(st_h5ad=..., cell_type_key=..., output_dir=..., plot_path=..., r_script_path=..., species=...)`
   - `species`: `"human"` or `"mouse"` — selects the CellChatDB; omission triggers gene-symbol detection.
   - Output: CellChat RDS object, communication/pathway tables, run settings and network plot.

## Outputs

- `results/08_neighborhood/<sample>_enrichment.csv` + heatmap PNG
- `results/09_communication/cellchat_network.png`
- `results/09_communication/cellchat_communication.csv` + `cellchat_pathways.csv` + `cellchat.rds` + `run_settings.json`

## Biological Interpretation

- Neighborhood: which cell-type pairs co-localize significantly (|z| > 2)? Does this match known tissue niches (e.g., tumor–macrophage, T cell–APC)?
- Communication: which ligand–receptor pairs dominate? Are they consistent with the tissue's known biology (e.g., immune checkpoint pairs in tumor)?

## Stop for Review

Present interpretation using the template from the parent `spatial-transcriptomics` skill. This is the final stage — after review, produce the overall workflow summary.
