---
name: spatial-visium-hd
description: Visium HD platform branch of the spatial transcriptomics workflow — reconstruct single cells from 2 μm bins via morphological segmentation and bin-to-cell aggregation. Use when the user's data has binned_outputs/square_002um (Visium HD Space Ranger output). Produces a cell-level h5ad, then stops for review.
license: MIT
---

# Visium HD Branch — Cell Segmentation & Bin-to-Cell

## Goal

Reconstruct **single-cell-level** expression from Visium HD 2 μm bins. The 2 μm bin is a transcript-capture unit, NOT a biological unit (a cell spans ~25 bins). The default 8 μm bin output is also NOT the analysis unit — **segment, not bin**.

## Prerequisites

- Visium HD Space Ranger output:
  - `binned_outputs/square_002um/filtered_feature_bc_matrix.h5` (bins × genes)
  - `binned_outputs/square_002um/spatial/tissue_positions.parquet` (bin coordinates)
  - `spatial/cytassist_image.tiff` (full-resolution HE image)
- Python: `scanpy`, `stardist`, `scipy`, `geopandas` (optional for soft assignment)

### GPU 加速（自动检测，无需手动配置）

- `segment_cells` 会自动探测 GPU：可用则用 GPU 加速，不可用自动回退 CPU。
- 有 NVIDIA GPU 时建议安装 GPU 版 TensorFlow：
  ```bash
  pip install "tensorflow[and-cuda]"
  ```
  无 GPU 或只想用 CPU：`pip install tensorflow` 即可，代码自动回退，不影响功能。
- 强制 CPU：设置环境变量 `STARDIST_DEVICE=cpu`。
- 注意：TF 首次在 GPU 上编译卷积会花 30-60 秒（XLA 编译），之后推理很快；大图分割 GPU 比 CPU 快约 6-10 倍。

## 参数决策约定（LLM 必须遵守）

**设计原则（2026-08-26）**：skill 不替 LLM 做"需要判断"的决策。参数分两类：

### A. 能直接计算的参数（代码自动算，LLM 不干预）
- `n_cells`（downsample 上限，每类型 1000）— 纯计算
- `n_neighbors`（表达图 KNN k）— 按数据规模 `min(15, max(5, n//1000))`
- `n_comps`（PCA 维数）— `min(30, n_cells-1, n_genes-1)`
- `HVG flavor`（counts vs log）— 数据自动检测（整数+非负=counts）
- `n_top_genes` — `min(2000, n_vars)`

### B. 需要 LLM 判断的参数（skill 输出诊断，LLM 决策后显式传参）
| 参数 | 为什么需要判断 | 判断依据 |
|---|---|---|
| `resolution`（域检测） | 依赖组织类型/分析粒度 | 小肠 10-20 域、脑皮层 7 层、肿瘤微环境可更多；skill 给出数据驱动默认 + 域数量诊断，LLM 看结果决定是否调整重跑 |
| `min_domains/max_domains` | 依赖组织先验 | 同上 |
| CellChat `trim` | 依赖数据稀疏度 | 检出率低用 0.1，高可收紧 |
| `expand_px`（核扩展） | 依赖组织细胞密度 | 默认 29px=8µm（10x 默认），LLM 看 overlay 判断 |
| 是否 destripe / Proseg | 依赖图像质量/用户需求 | 看 03_nuclei_crop.png 判断 |

**工作流**：
1. skill 先跑可计算参数（A 类）
2. 需要判断处（B 类）→ skill 输出 DIAGNOSTIC 信息（域数量、数据规模等）
3. **LLM 看到诊断后，根据组织学知识/用户需求决策**，显式传参重跑（如 `resolution=1.2`）
4. 禁止 skill 内部用参数搜索"自动找正确答案"——那是不可判定的

## Steps

流程按 stage 编号执行（S1→S6），S4 分两个分支：

1. **S1 Load 2 μm bins**
   - `sc.read_10x_h5(square_002um/filtered_feature_bc_matrix.h5)` — millions of bins, keep sparse.
   - Attach coordinates from `tissue_positions.parquet` into `obsm['spatial']`.

2. **S2 Destripe (recommended)**
   - Correct the known Visium HD striping artifact (per-row/column quantile scaling of UMI counts at 2 μm).
   - Optional but improves segmentation downstream.

3. **S3 Segment nuclei (StarDist)**
   - Model: `2D_versatile_he` (H&E) on the full-res HE image.
   - Tiled inference (`predict_instances_big`, block_size ~4096, overlap ~128).
   - Save labels as sparse NPZ. **Inspect a crop** (`03_nuclei_crop.png`) to confirm segmentation quality before proceeding.

4. **S4 Cell reconstruction（两分支二选一或都跑）**
   - **S4A StarDist + bin2cell（默认）**: 核多边形 buffer 扩展（`expand_px` 默认 29px ≈ 8µm，对齐 10x Space Ranger `--nucleus-expansion-distance-micron` 默认值；0 = 严格核内）→ bin-to-cell 聚合（重叠 bins 按最近核分配）→ `04_bin2cell.h5ad`。
   - **S4B Proseg**: run Proseg Bayesian segmentation (voxel_size auto, samples auto, 2 μm bins + StarDist masks) → `proseg_out.zarr/` + `04_proseg_cells.h5ad`。
   - 两分支共用 S3 的核分割结果；输出后各生成一张全图高清几何 overlay（`04_cells_overlay.png`）。
   - **分割完成后得到的是"空间单细胞"数据，直接按 scRNA-seq 处理（注释 → domains → neighborhood → CellChat），不需要反卷积**（依据：Bin2cell/ENACT/10x 官方指南，见知识库 `concepts/visium-hd-segmentation-vs-deconvolution.md`）。
   - 细胞类型鉴定用**注释（annotation / label transfer）**：CellTypist、CellAssign、或 scRNA 参考（如 Haber 2017 小鼠小肠）做 label transfer，不用 deconvolve。
   - **反卷积（deconvolve）不是 10x 官方流程的一部分**（无论分割与否）：官方做法是聚类 + marker 注释。反卷积只是社区工具（Cell2location/SpatialDWLS/SpaceXR）的可选增强，且需要外部 scRNA 参考——默认不做，用户明确要求时才走经典 visium.py 的 `deconvolve_spatial_*`。

5. **S5 QC & sanity check**
   - Report: number of cells, median genes/cell, median counts/cell, median bins/cell, novelty score, fraction of tissue bins assigned.
   - Flag: suspiciously low cell count, high empty fraction, low bins/cell (median < 5 → 提示 bin_to_cell 缺核扩展), or segmentation failures.
   - 输出 `segmentation_summary.json` + `qc_metrics.json` + **QC 可视化面板**（`05_qc.png`：基因数/UMI 数/mito% violin + novelty score/bin_count 直方图 + 空间低质量细胞标记，参照 SIB-Swiss 空间转录组培训与 bcbio spatial-reports 标准）+ 全图 overlay；**提示用户可用 ROI 工具框选局部检查**（见下方"HD 可视化 & 交互式 ROI 选择"）。

## Outputs

统一输出到 `results/07b_segmentation/`，每个方法一个子目录：

```
results/07b_segmentation/
├── <sample>_stardist/            # StarDist 分支
│   ├── 01_load_bins.h5ad
│   ├── 02_destripe.h5ad          # (如启用 destripe)
│   ├── 03_nuclei_labels.npz      # 核分割 labels
│   ├── 03_nuclei_polys.pckl      # StarDist 多边形 (供几何绘图)
│   ├── 03_nuclei_crop.png        # 核分割检查图
│   ├── 04_bin2cell.h5ad          # 细胞级 AnnData (cells × genes)
│   ├── 04_cells_overlay.png      # 全图高清几何 overlay (蓝色轮廓)
│   ├── 05_qc.png                 # QC 可视化面板 (基因/UMI/bin 分布 + 空间)
│   ├── segmentation_summary.json # 细胞数/QC/参数
│   ├── qc_metrics.json           # QC 数字指标 (n_cells/median genes/counts)
│   └── params.json               # 实际运行参数
├── <sample>_proseg/              # Proseg 分支
│   ├── 01_load_bins.h5ad
│   ├── 03_nuclei_labels.npz / 03_nuclei_polys.pckl
│   ├── 03_nuclei_crop.png
│   ├── proseg_out.zarr/          # Proseg 几何输出 (cell_boundaries)
│   ├── 04_proseg_cells.h5ad      # 细胞级 AnnData
│   ├── 04_cells_overlay.png      # 全图高清几何 overlay (红色轮廓)
│   ├── segmentation_summary.json
│   └── params.json
└── crops/                        # 用户 ROI 选择产出的 crop 图
    ├── roi_<sample>_<method>.json  # ROI 坐标 (pick_roi.py 输出)
    └── crop_<sample>_<method>_<x0>_<y0>_<x1>_<y1>.png
```

统一命名规则：
- 中间产物：`NN_<阶段>.<ext>`（01_load / 02_destripe / 03_nuclei / 04_cells）
- 图：`04_cells_overlay.png`（全图高清几何）、`03_nuclei_crop.png`（核分割检查）
- 分支标识：目录名 `<sample>_stardist` / `<sample>_proseg`
- 所有用户 ROI 产出集中在 `crops/` 子目录

## 下游验证（S4 之后、交用户审查前）

细胞级 h5ad 产出后，先跑下游验证确认结果可进入共享下游流程，再交用户审查：

```bash
# 服务器 (用 S4 输出的细胞级 h5ad)
python validate_downstream.py \
    --input <sample>_stardist/04_bin2cell.h5ad \
    --out-dir downstream_test \
    --n-cells 10000 --n-top-genes 2000 --n-pcs 30 --resolution 1.0
```

输出（默认 `downstream_test/` 目录）：
- `sub10000.h5ad` — 验证子集（原始计数）
- `sub10000_normalized.h5ad` + `sub10000_hvg.png` — normalize + 高变基因图
- `sub10000_clustered.h5ad` + `sub10000_umap.png` — PCA + Leiden 聚类 + UMAP 图

全部输出存在才算通过；任一项缺失脚本报错退出。验证通过后，再向用户展示全图 overlay 并提示可用 ROI 工具做局部检查。

## HD 可视化 & 交互式 ROI 选择

绘图统一走 `visium_new_platforms.py` 中的 `_plot_cells_overlay` / `plot_cells_crop`，输出 **HE 全分辨率 + 真实细胞几何（多边形边界）** 的高清图，不再是低清质心散点。

### 挂载点：S4 分析结束后

**S4 输出结果和大图后，提示用户：可以用 ROI 工具选择区域做局部检查。**

```text
[Visium HD] S4 完成:
  StarDist: N cells -> <sample>_stardist/04_cells_overlay.png (全图蓝色轮廓)
  Proseg:   M cells -> <sample>_proseg/04_cells_overlay.png   (全图红色轮廓)
  如对某区域细胞边界与 HE 组织对齐存疑，可用 ROI 工具框选局部放大检查:
    python pick_roi.py --image tissue_image.png --out crops/roi_<sample>_<method>.json
    python make_crop.py --method stardist|proseg --roi-json crops/roi_<sample>_<method>.json
```

### 全图高清几何 overlay

- StarDist 分支：蓝色轮廓（`edge_color="tab:blue"`，linewidth 0.12），输出 `<sample>_stardist/04_cells_overlay.png`
- Proseg 分支：红色轮廓（`edge_color="tab:red"`，linewidth 0.12），输出 `<sample>_proseg/04_cells_overlay.png`
- 默认：HE 全分辨率（downsample=1）、全部细胞（max_points=None）、dpi=150

### 交互式 ROI 选择（用户手动框选区域）

本地工具 `pick_roi.py`（弹窗 + 鼠标框选），输出到 `results/07b_segmentation/crops/`：

```bash
# 本地 (需 Python + matplotlib + pillow)
python pick_roi.py --image tissue_image.png --out roi_<sample>_<method>.json
# 鼠标左键拖拽框选 ROI; Enter 确认, r 重选, q 退出
# 输出全分辨率像素坐标到 roi_<sample>_<method>.json
```

服务器配套 `make_crop.py`（用 ROI 坐标生成高清 crop，输出到 `results/07b_segmentation/crops/`）：

```bash
# 服务器 (crop 图统一存 crops/ 子目录)
python make_crop.py --method stardist --roi x0,y0,x1,y1 [--out crops/crop_<sample>_stardist_x0_y0_x1_y1.png]
python make_crop.py --method proseg  --roi-json crops/roi_<sample>_proseg.json
```

`plot_cells_crop` 函数本身支持三种 ROI 指定方式：
- `roi=None`：自动选细胞最密集区
- `roi=(x0, y0, x1, y1)`：手动矩形区域
- `roi=(cx, cy)` + `roi_size`：中心点 + 边长

## 空间域命名规范（LLM 必须遵守）

**背景（2026-08-26）**：空间转录组领域没有统一的 domain 命名标准（不像细胞类型有 Cell Ontology / CellTypist 参考库）。命名按以下 5 步流程执行，**禁止 LLM 自由发挥**。

### LLM 命名流程（5 步）

1. **读 annotation CSV**：`identify_spatial_domains` 产出的 `*_domains_annotation.csv`，含每 domain 的 top marker 基因（logFC）+ 细胞类型组成（占比）
2. **判定组织类型**：从样本元数据/组织来源推断（肠/肺/心/肝/癌等）。不同组织的解剖结构差异大，命名词汇必须匹配组织
3. **对照组织类型参考表选候选名**：优先解剖结构，其次组织学特征。参考表（持续积累，遇到新组织补充）：
   - 小肠：绒毛吸收上皮 / 隐窝-绒毛过渡 / 隐窝底部(潘氏细胞) / 黏膜下间质 / 平滑肌 / 淋巴组织(HEV+) / B细胞区 / 浆细胞富集区
   - 脑：皮层 L1-L6 / 白质 / 海马区 / 脑室周围（参照 DLPFC spatialLIBD 标准）
   - 肺：肺泡区 / 支气管上皮 / 血管周围间质 / 平滑肌 / 淋巴滤泡 / 肿瘤实性区 / 坏死区
   - 心脏：心肌层 / 心外膜 / 心内膜 / 血管壁 / 纤维化区 / 脂肪浸润区
   - 肝脏：肝小叶中央 / 门静脉周围 / 胆管区 / 纤维化区 / 肿瘤结节
   - 癌症（通用）：肿瘤实性区 / 肿瘤浸润前沿 / 肿瘤间质 / 免疫浸润区 / 坏死区 / 淋巴聚集区 / 血管区 / 纤维包膜
4. **用 marker 基因验证**：候选名必须有 marker 支持。例："隐窝底部"需 Lgr5/Defa 等潘氏/干细胞 marker；"心肌层"需 TNNT2/ACTA2 等；"肿瘤浸润前沿"需 EMT/增殖 marker
5. **输出统一格式**：`D{n} {结构名}({特征})`。证据不足 → `D{n} 未命名(待定)` 并说明缺什么证据（如"无明确 marker 支持，需 H&E 图像辅助判断"）

### 命名约束

- **domain 是组织切片上的空间结构域，不是细胞类型**——命名必须回答"这个区域在组织里是什么"（解剖结构/组织学特征），细胞类型组成只作为佐证
- 命名对齐 Uberon 解剖学术语（跨物种解剖学本体），不自己造词
- 命名后向用户展示：annotation CSV 摘要 + 拟定命名 + 依据（marker/细胞类型），**用户确认后才画图**
- 用户指出命名错误 → 记录到 `knowledge/lessons/`，并补充组织类型参考表

## Biological Interpretation

- Report total cells and QC stats.
- Verify reconstructed cells overlap tissue regions in the HE overlay; flag low-density regions.
- Cross-check: do cell-level patterns respect tissue architecture (e.g., epithelial sheets, immune infiltrates)?

## Stop for Review

Present interpretation using the template from the parent `spatial-transcriptomics` skill. Wait for `通过` / `调整` / `跳过` before proceeding to shared downstream.

## Notes

- Memory: 2 μm matrices are huge (millions of bins) — prefer sparse/chunked operations.
- Deconvolution is NOT needed for Visium HD — cells are already resolved after bin-to-cell.
- Reference implementations: bin2cell (Teichmann lab), ENACT (Sanofi) — in `knowledge/references/projects/spatial-transcriptomics-seg/`.
