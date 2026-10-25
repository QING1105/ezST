---
name: spatial-visium
description: Classic Visium platform branch of the spatial transcriptomics workflow — spot-level 5-stage pipeline (load+QC, normalize+cluster, domains+SVG, deconvolution, downstream prep). Use when the user's data is classic 10x Visium spot-level output (55 μm spots). Produces a spot-level h5ad + deconvolution results, then stops for review.
license: MIT
---

## 0.2.1 运行约定

表达聚类使用 `identify_spatial_domains(..., spatial_weight=0)`；加入空间邻接时显式设置 0 到 1 之间的权重，并报告 `uns['domain_clustering']`。CellChat 默认 `mode="expression"`；空间模式需要 `mode="spatial"`、`microns_per_coordinate`、`spatial_tolerance`，按样本尺度设置作用范围和距离缩放。输入优先使用保存的原始 counts。CellChat 输出实际为 `cellchat_communication.csv`、`cellchat_pathways.csv`、`cellchat.rds`、`cellchat_network.png` 与 `run_settings.json`。

StarDist 几何使用数值 `.npz`，旧 pickle 读取需显式 `allow_legacy_pickle=True`。HD 裁图通过 `--image`、`--cells`、`--geometry`、`--roi` 指定输入；全部坐标使用全分辨率图像像素。Proseg 基因名来自其导出的 `proseg_genes.csv.gz`，导入后核对 `var_names`。Python 与插件内运行文件通过 `python scripts/sync_plugin_runtime.py` 同步。



# Classic Visium Branch — Spot-Level 5-Stage Pipeline

## Goal

Run the classic Visium spot-level analysis. Spots are 55 μm and cover multiple cells — **deconvolution IS required** to estimate cell-type composition (unlike Visium HD / Xenium / Atera which are cell-resolved).

## Stages

### S1 — Load + QC
- Tool chain: `init_spatial_project` → `load_visium_data` → `filter_visium_spots`
- Plots: QC violin plots / spot count distributions before and after filtering
- Results: filtered h5ad, QC summary
- Review focus: are the filtering thresholds (`min_counts`, `min_genes`, `pct_mt`) reasonable? Is the retained spot count consistent with the expected tissue size?

### S2 — Normalize + Cluster
- Tool chain: `normalize_visium` → `cluster_spatial_data`
- Plots: UMAP, spatial scatter colored by cluster
- Results: normalized/clustered h5ad
- Review focus: is the number of clusters biologically plausible? Do clusters segregate spatially?

### S3 — Spatial Domains + SVG
- Tool chain: `identify_spatial_domains` → `find_spatially_variable_genes`
- Plots: spatial domain map, top SVG spatial plots
- Results: domain assignments, SVG ranked CSV
- Review focus: do spatial domains match known tissue architecture? Do top SVGs make biological sense?

### S4 — Deconvolution (required for classic Visium)
- Tool chain: `deconvolve_spatial_destvi` (preferred) or `deconvolve_spatial_spotlight`
- Requires reference scRNA-seq h5ad with cell-type labels in `obs`
- Plots: cell-type proportion spatial maps, proportion stacked bar
- Results: proportions CSV per spot, deconvolved h5ad
- Review focus: are the dominant cell types consistent with the tissue type? Any unexpected cell type dominating?

### S5 — Prepare for Shared Downstream
- Ensure spot-level h5ad has cell-type labels (from deconvolution) attached for neighborhood enrichment / CellChat.
- If the user also wants single-cell reconstruction from classic Visium (optional, unusual), mention STIE as the morphological-deconvolution route — but this is NOT the default.

## Outputs

- `results/01_loading/` ... `results/07_deconvolution/` (see parent skill conventions)
- Final: spot-level h5ad with deconvolution proportions attached

## Biological Interpretation

Use the parent skill's interpretation template for each stage. For S4 specifically, cross-check dominant cell types against tissue type (e.g., gastric cancer should show epithelial/tumor cells in tumor regions).

## Stop for Review

Present interpretation at every stage using the template from the parent `spatial-transcriptomics` skill. Wait for `通过` / `调整` / `跳过` before proceeding.
