---
name: spatial-transcriptomics
description: End-to-end 10x Visium spatial transcriptomics analysis workflow with staged execution and human review gates. Use when the user asks to analyze Visium / spatial transcriptomics data, run spatial analysis pipelines (load, QC, normalize, cluster, spatial domains, SVG, deconvolution, neighborhood enrichment, cell-cell communication), or interpret spatial plots biologically. Each stage produces plots and results, the LLM interprets their biological meaning, then stops and waits for the user's review before proceeding.
license: MIT
---

# Spatial Transcriptomics (10x Visium) Staged Workflow

## Overview

Run 10x Visium spatial transcriptomics analysis as a **staged pipeline**. The workflow is split into 5 stages. Each stage:

1. **Run** the tool chain for that stage, producing plots (PNG) and result files.
2. **Interpret** the outputs with the vision tool (for plots) and by reading result files. Write a biological interpretation using the template below.
3. **Stop** and present the interpretation + plots to the user.
4. **Wait for review**. The user chooses one of: `通过` (proceed), `调整` (adjust parameters and rerun this stage), or `跳过` (skip this stage and continue).

Do NOT proceed to the next stage without an explicit user decision.

## Environment Check (run first)

Before any analysis, verify the environment is ready. If not, guide the user to install it (one command, no manual steps):

```bash
# Check critical Python packages
python -c "import scanpy, squidpy, scvi, leidenalg; print('OK')"
# Check critical R packages (only needed for S4/S5)
Rscript -e "library(Seurat); library(SPOTlight); library(CellChat); cat('OK\n')"
```

If any Python import fails:
1. Tell the user ezST is not installed.
2. Install it with pip (one command, auto-installs all Python deps):
   ```bash
   pip install git+https://github.com/QING1105/ezST.git
   ```
3. Re-run the Python check above.

If R packages are missing (only needed for S4 deconvolution / S5 communication):
1. Tell the user the R deps are missing.
2. Run the R installer (single command):
   ```bash
   Rscript install/install_r_packages.R
   ```
   (needs R >= 4.3; Windows needs Rtools)
3. Re-run the R check above.

Do NOT start analysis until the required checks pass (S1-S3 only need Python; S4/S5 also need R).

## Data Requirements

The workflow accepts three input formats (auto-detected):

| Format | Files | Notes |
|:--|:--|:--|
| h5ad | single `.h5ad` with `obsm['spatial']` | Preferred; coordinates embedded |
| 10x matrix | `matrix.mtx.gz` + `barcodes.tsv.gz` + `features.tsv.gz` | Standard CellRanger output |
| wide CSV | counts CSV + spot coordinates CSV | Rows = barcodes, cols = genes |

Reference scRNA-seq (h5ad with cell-type labels in `obs`) is required for deconvolution stages.

## Stage Protocol

### S1 — Load + QC
- Tool chain: `init_spatial_project` → `load_visium_data` → `filter_visium_spots`
- Plots: QC violin plots / spot count distributions before and after filtering
- Results: filtered h5ad, QC summary
- Review focus: are the filtering thresholds (`min_counts`, `min_genes`, `pct_mt`) reasonable? Is the retained spot count consistent with the expected tissue size?

### S2 — Normalize + Cluster
- Tool chain: `normalize_visium` → `cluster_spatial_data`
- Plots: UMAP, spatial scatter colored by cluster
- Results: normalized/clustered h5ad
- Review focus: is the number of clusters biologically plausible? Do clusters segregate spatially (not just noise)?

### S3 — Spatial Domains + SVG
- Tool chain: `identify_spatial_domains` → `find_spatially_variable_genes`
- Plots: spatial domain map, top SVG spatial plots
- Results: domain assignments, SVG ranked CSV
- Review focus: do spatial domains match known tissue architecture (e.g., tumor regions, immune infiltrate)? Do top SVGs make biological sense for the tissue?

### S4 — Deconvolution
- Tool chain: `deconvolve_spatial_destvi` (preferred) or `deconvolve_spatial_spotlight`
- Plots: cell-type proportion spatial maps, proportion stacked bar
- Results: proportions CSV per spot, deconvolved h5ad
- Review focus: are the dominant cell types consistent with the tissue type? Any unexpected cell type dominating?

### S5 — Downstream (Neighborhood + Communication)
- Tool chain: `spatial_neighborhood_enrichment` → `infer_spatial_cell_communication`
- Plots: neighborhood enrichment heatmap, CellChat network/heatmap plots
- Results: enrichment z-score table, CellChat object + summary
- Review focus: are enriched co-localizations and inferred ligand-receptor interactions biologically plausible?

## Biological Interpretation Template

Use this exact structure for every stage's interpretation:

```
【S{n} 阶段结果】
- 关键数字：XXX spots retained / XXX clusters / XXX spatial domains / top cell types...
- 图片观察：图中可见……（描述 spatial pattern, e.g., 高表达区域聚集在左上象限）
- 生物学意义：这提示……（e.g., 该区域可能是肿瘤浸润边缘，富集免疫细胞通讯）
- 建议：是否调整参数 / 进入下一阶段
```

## Review Protocol

After presenting the interpretation, ask the user explicitly (do not assume):

- `通过` → proceed to next stage
- `调整` → ask which parameters to change, rerun only the current stage
- `跳过` → skip current stage, continue to next

Record the user's decision for each stage in the final summary.

## Output Conventions

- All results under the project directory created by `init_spatial_project` (e.g., `/data/gastric/proj/`).
- Naming: `results/01_loading/`, `results/02_qc/`, `results/03_normalization/`, `results/04_clustering/`, `results/05_domains/`, `results/06_svg/`, `results/07_deconvolution/`, `results/08_neighborhood/`, `results/09_communication/`.
- Final summary at the end: one line per stage with user decision (通过/调整/跳过).
