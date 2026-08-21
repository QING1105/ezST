# Spatial Transcriptomics (10x Visium)

Staged 10x Visium spatial transcriptomics analysis for Codex, with human review gates at every stage.

## Skills included

- **spatial-transcriptomics** — orchestrator: 5-stage protocol + biological interpretation template + review protocol (通过 / 调整 / 跳过)
- **spatial-s1-load-qc** — load Visium data, QC-filter spots
- **spatial-s2-normalize-cluster** — normalize, log-transform, cluster (PCA + UMAP + Leiden)
- **spatial-s3-domains-svg** — spatial domains, Moran's I SVG detection
- **spatial-s4-deconvolve** — DestVI / SPOTlight cell-type deconvolution
- **spatial-s5-downstream** — neighborhood enrichment + CellChat communication

## Example prompts

```text
Analyze my gastric cancer Visium data at '/data/gastric/GSE246011/' with the staged workflow, stopping at each stage for my review.
```

```text
Deconvolve my Visium spatial data at '/data/gastric/proj/results/01_loading/sample.h5ad' using the scRNA-seq reference at '/data/ref_scRNA.h5ad'.
```

```text
Run cell-cell communication analysis on my clustered Visium data at '/data/gastric/proj/results/04_clustering/sample_clustered.h5ad'.
```
