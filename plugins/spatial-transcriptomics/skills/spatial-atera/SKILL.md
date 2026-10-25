---
name: spatial-atera
description: Atera platform branch of the spatial transcriptomics workflow — load and validate Atera cell-level output (AnnData + Zarr segmentation) for downstream analysis. Use when the user's data is from the 10x Atera In Situ instrument (whole-transcriptome in situ, 2026). Produces a cell-level h5ad, then stops for review.
license: MIT
---

## 0.2.1 运行约定

表达聚类使用 `identify_spatial_domains(..., spatial_weight=0)`；加入空间邻接时显式设置 0 到 1 之间的权重，并报告 `uns['domain_clustering']`。CellChat 默认 `mode="expression"`；空间模式需要 `mode="spatial"`、`microns_per_coordinate`、`spatial_tolerance`，按样本尺度设置作用范围和距离缩放。输入优先使用保存的原始 counts。CellChat 输出实际为 `cellchat_communication.csv`、`cellchat_pathways.csv`、`cellchat.rds`、`cellchat_network.png` 与 `run_settings.json`。

StarDist 几何使用数值 `.npz`，旧 pickle 读取需显式 `allow_legacy_pickle=True`。HD 裁图通过 `--image`、`--cells`、`--geometry`、`--roi` 指定输入；全部坐标使用全分辨率图像像素。Proseg 基因名来自其导出的 `proseg_genes.csv.gz`，导入后核对 `var_names`。Python 与插件内运行文件通过 `python scripts/sync_plugin_runtime.py` 同步。



# Atera Branch — Cell-Level Loading & Validation

## Goal

Load Atera data as a **cell-level** matrix. Atera is 10x's next-generation in situ platform (launched 2026) that outputs whole-transcriptome spatial data at single-cell sensitivity with **platform-provided segmentation** — no re-segmentation and no deconvolution needed.

## Prerequisites

- Atera output (official formats):
  - Cell-feature matrix: **AnnData** (cells × genes)
  - Transcripts and segmentation boundaries: **Zarr**
  - Tissue images: **OME-TIFF**
- Python: `scanpy`, `zarr` (if reading segmentation)

## Steps

1. **Load cell matrix**
   - AnnData: `ad.read_h5ad(...)` — cells × genes.
   - If only Zarr + transcripts provided, aggregate transcripts to cells using segmentation boundaries first (rare; usually AnnData is provided).

2. **Validate platform segmentation**
   - Confirm cell boundaries exist (from Atera pipeline).
   - Do NOT re-segment unless the user explicitly asks.

3. **QC & sanity check**
   - Report: number of cells, median genes/cell, median counts/cell.
   - Flag: very low counts, very high counts, abnormal spatial distribution.

## Outputs

- `results/01_loading/<sample>_cells.h5ad` — cell-level AnnData
- `results/01_loading/QC_plots.png` — QC violin plots

## Biological Interpretation

- Report total cells and QC stats.
- Note tissue type and expected cell composition.
- Atera is whole-transcriptome — check gene coverage (18,000+ expected).

## Stop for Review

Present interpretation using the template from the parent `spatial-transcriptomics` skill. Wait for `通过` / `调整` / `跳过` before proceeding to shared downstream.

## Notes

- Atera segmentation comes with the platform — treat as ground truth unless the user asks otherwise.
- Deconvolution is NOT needed for Atera — cells are already resolved.
- Atera is new (2026); formats may evolve — check the instrument's current output spec if files differ from the above.
