# Spatial Transcriptomics (10x Visium / Visium HD)

Staged spatial transcriptomics analysis for Codex, with human review gates at every stage.

## Supported platforms & branches

| Platform | Skill | Analysis unit | Branch |
|:--|:--|:--|:--|
| 10x Visium (v1/v2) | `spatial-transcriptomics` (classic) | spot | load+QC → normalize+cluster → domains+SVG → **deconvolve** → neighborhood+CellChat |
| **10x Visium HD** | `spatial-visium-hd` | **single cell** (after segmentation) | segment (bin2cell / Proseg) → QC → **annotate (label transfer)** → domains → neighborhood+CellChat |

> Verified platforms: Classic Visium and Visium HD. Other platforms are NOT yet validated in this repo.
> Visium HD note: after segmentation the data is single-cell level, so it is processed **as scRNA-seq** (annotation via label transfer), **not deconvolved**. See `spatial-visium-hd/SKILL.md`.

## Skills included

- **spatial-transcriptomics** — orchestrator: 5-stage protocol + biological interpretation template + review protocol (通过 / 调整 / 跳过)
- **spatial-s1-load-qc** — load Visium data, QC-filter spots
- **spatial-s2-normalize-cluster** — normalize, log-transform, cluster (PCA + UMAP + Leiden)
- **spatial-s3-domains-svg** — spatial domains, Moran's I SVG detection
- **spatial-s4-deconvolve** — DestVI / SPOTlight cell-type deconvolution
- **spatial-s5-downstream** — neighborhood enrichment + CellChat communication
- **spatial-visium-hd** — Visium HD branch: bin2cell / Proseg segmentation → single-cell downstream

## Changelog

- **0823**: fixed a bug in the uploaded content; re-pushed to git.
- **0826**: added Visium HD support — two segmentation methods (10x official bin2cell and Proseg); downstream refactored to LLM-in-the-loop decisions (auto species detection, data-driven parameters). See the full pipeline figure in the main README (iterations, 0826).

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

```text
Analyze my Visium HD sample at '/data/visium_hd/' — segment cells (bin2cell), annotate with Haber 2017 reference, then run domains and CellChat.
```

