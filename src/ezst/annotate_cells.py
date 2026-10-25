#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""annotate_cells.py — Visium HD 分割后细胞注释（双方法）

支持两种 label transfer 方法，按参考数据格式自动选择或手动指定：

  A. scVI（深度生成模型）: 要求参考数据是**原始 counts**（整数、无负值）。
     Haber 原始版（GEO 下载自建 h5ad）等适用。
  B. PCA + kNN ingest（Seurat 式）: 不要求原始 counts，直接处理 log1p/归一化数据。
     Haber Besca 处理版（Zenodo 4447233，log1p）等适用。

方法选择（--method auto，默认）:
  - 检测参考数据 X: 全整数且 min>=0 -> scVI
  - 否则（有负值/非整数/log1p）-> ingest

输入:
  --cells    Visium HD 细胞级 h5ad (bin2cell/proseg 输出, 有 obsm['spatial'])
  --ref      参考 scRNA h5ad (有 obs[cell_type_key])
  --cell-type-key  参考中的细胞类型列
  --method   auto | scvi | ingest (默认 auto)
  --out-dir  输出目录

流程:
  1. 基因取交集, 统一数据形态 (scVI: 需 counts; ingest: 需 log1p)
  2. scVI: 拼接 -> scVI 训练 -> 潜空间 kNN 预测
     ingest: 拼接 -> HVG -> PCA -> sklearn kNN 预测
  3. 输出: 带注释的细胞 h5ad + 注释统计图 + 空间分布图

依据: 分割后得到"空间单细胞", 用 label transfer 注释, 不走 deconvolve
  (见 knowledge/concepts/visium-hd-segmentation-vs-deconvolution.md)
"""
import argparse
import os
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
import scipy.sparse as sp


def _is_counts(X):
    """判断矩阵是否原始 counts: 稀疏取 data, 检查全整数且 min>=0"""
    if sp.issparse(X):
        d = X.data
    else:
        d = np.asarray(X).ravel()
    if len(d) == 0:
        return True
    return bool(d.min() >= 0 and np.allclose(d, np.round(d)))


def _to_log1p(a):
    """若是 counts, 转 log1p (ingest 分支用)"""
    if not _is_counts(a.X):
        return a
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def _to_counts(a):
    """Use preserved raw counts for the count-likelihood model."""
    try:
        from .visium import _raw_count_adata
    except ImportError:
        from visium import _raw_count_adata
    return _raw_count_adata(a, "annotation input")


def run_scvi(cells_s, ref_s, cell_type_key, out_dir, args):
    """scVI label transfer (需原始 counts)"""
    import scvi
    print("[scVI] training...")
    cells_s.obs["_batch"] = "query"
    ref_s.obs["_batch"] = "ref"
    combined = ad.concat([ref_s, cells_s], join="outer", merge="same")
    combined = combined[:, cells_s.var_names.intersection(ref_s.var_names)].copy()
    combined.obs["_batch"] = combined.obs["_batch"].astype(str)
    scvi.model.SCVI.setup_anndata(combined, batch_key="_batch")
    model = scvi.model.SCVI(combined, n_latent=30, n_layers=2)
    model.train(max_epochs=args.max_epochs)
    combined.obsm["X_scVI"] = model.get_latent_representation()
    # kNN on latent space
    from sklearn.neighbors import KNeighborsClassifier
    ref_mask = combined.obs["_batch"] == "ref"
    query_mask = ~ref_mask
    knn = KNeighborsClassifier(n_neighbors=min(10, int(ref_mask.sum())), weights="distance", metric="cosine")
    knn.fit(combined.obsm["X_scVI"][ref_mask], combined.obs.loc[ref_mask, cell_type_key].astype(str).to_numpy())
    pred = knn.predict(combined.obsm["X_scVI"][query_mask])
    proba = knn.predict_proba(combined.obsm["X_scVI"][query_mask])
    return pred, proba


def run_ingest(cells_s, ref_s, cell_type_key, out_dir, args):
    """PCA + kNN label transfer (Seurat 式, 不需 counts)"""
    from sklearn.neighbors import KNeighborsClassifier
    print("[ingest] HVG + PCA...")
    cells_s.obs["_batch"] = "query"
    ref_s.obs["_batch"] = "ref"
    combined = ad.concat([ref_s, cells_s], join="outer", merge="same")
    common = cells_s.var_names.intersection(ref_s.var_names)
    combined = combined[:, common].copy()
    combined.obs["_batch"] = combined.obs["_batch"].astype(str)
    try:
        sc.pp.highly_variable_genes(combined, n_top_genes=min(args.n_top_genes, len(common) - 1))
    except Exception:
        combined.var["highly_variable"] = True
    combined = combined[:, combined.var["highly_variable"]].copy()
    sc.pp.scale(combined, max_value=10)
    sc.tl.pca(combined, n_comps=min(30, combined.n_obs - 1, combined.n_vars - 1), svd_solver="arpack")
    print("[ingest] kNN label transfer...")
    ref_mask = combined.obs["_batch"] == "ref"
    query_mask = ~ref_mask
    knn = KNeighborsClassifier(n_neighbors=min(10, int(ref_mask.sum())), weights="distance", metric="cosine")
    knn.fit(combined.obsm["X_pca"][ref_mask], combined.obs.loc[ref_mask, cell_type_key].astype(str).to_numpy())
    pred = knn.predict(combined.obsm["X_pca"][query_mask])
    proba = knn.predict_proba(combined.obsm["X_pca"][query_mask])
    return pred, proba


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--cell-type-key", default="cell_type")
    ap.add_argument("--method", default="auto", choices=["auto", "scvi", "ingest"])
    ap.add_argument("--out-dir", default="./annotation_results")
    ap.add_argument("--n-top-genes", type=int, default=2000)
    ap.add_argument("--max-epochs", type=int, default=80)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")

    # --- 1. 加载 ---
    print("Loading cells...")
    cells = ad.read_h5ad(args.cells)
    print("Loading ref...")
    ref = ad.read_h5ad(args.ref)
    if args.cell_type_key not in ref.obs.columns:
        for cand in ["cell_type", "cell_ontology_class", "annotation", "labels", "celltype", "celltype1"]:
            if cand in ref.obs.columns:
                args.cell_type_key = cand
                break
    print(f"ref cell_type_key: {args.cell_type_key}")
    print(f"ref types: {ref.obs[args.cell_type_key].nunique()}, cells: {ref.n_obs}")
    print(f"cells: {cells.n_obs} x {cells.n_vars}")

    # --- 2. 数据格式诊断 + 方法选择 ---
    ref_is_counts = _is_counts(ref.X)
    cells_is_counts = _is_counts(cells.X)
    print(f"ref X is counts: {ref_is_counts} (min={ref.X.min() if not sp.issparse(ref.X) else ref.X.data.min():.2f}, max={ref.X.max() if not sp.issparse(ref.X) else ref.X.data.max():.2f})")
    print(f"cells X is counts: {cells_is_counts}")

    if args.method == "auto":
        method = "scvi" if ref_is_counts else "ingest"
    else:
        method = args.method
    print(f"method: {method}")

    # --- 3. 基因交集 ---
    common = cells.var_names.intersection(ref.var_names)
    print(f"common genes: {len(common)}")
    ref_s = ref[:, common].copy()
    cells_s = cells[:, common].copy()

    # --- 4. 数据形态统一 ---
    if method == "scvi":
        # scVI 需要 counts
        print("scVI: loading preserved raw counts")
        cells_s = _to_counts(cells_s)
        ref_s = _to_counts(ref_s)
        pred, proba = run_scvi(cells_s, ref_s, args.cell_type_key, args.out_dir, args)
    else:
        # ingest 需要 log1p
        print("ingest: converting to log1p (if needed)")
        cells_s = _to_log1p(cells_s)
        ref_s = _to_log1p(ref_s)
        pred, proba = run_ingest(cells_s, ref_s, args.cell_type_key, args.out_dir, args)

    # --- 5. 输出 ---
    cells.obs["cell_type"] = pred
    cells.obs["cell_type_prob"] = proba.max(axis=1)
    out_h5ad = os.path.join(args.out_dir, "annotated_cells.h5ad")
    cells.write(out_h5ad, compression="lzf")
    print(f"saved -> {out_h5ad}")

    # 注释统计图
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    counts_ = cells.obs["cell_type"].value_counts()
    axes[0].barh(counts_.index, counts_.values, color="#4C72B0")
    axes[0].set_title(f"Cell type counts ({method}, n={cells.n_obs:,})")
    axes[0].set_xlabel("cells")
    axes[1].hist(cells.obs["cell_type_prob"], bins=50, color="#55A868")
    axes[1].set_title("Annotation confidence (kNN proba)")
    axes[1].set_xlabel("probability")
    plt.tight_layout()
    png = os.path.join(args.out_dir, "annotation_summary.png")
    plt.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # 空间分布图
    if "spatial" in cells.obsm:
        coords = cells.obsm["spatial"]
        fig, ax = plt.subplots(1, 1, figsize=(12, 12))
        types = cells.obs["cell_type"].astype("category")
        # FIX 2026-08-26: tab20 只有 20 色, >20 类会重复. 用统一 _categorical_cmap.
        from visium import _categorical_cmap
        cmap_colors = _categorical_cmap(types.cat.categories.size)
        colors_map = {c: cmap_colors[i] for i, c in enumerate(types.cat.categories)}
        if len(coords) > 100000:
            rng = np.random.default_rng(0)
            idx = rng.choice(len(coords), 100000, replace=False)
            coords_s, types_s = coords[idx], types.iloc[idx]
        else:
            coords_s, types_s = coords, types
        for c in types_s.cat.categories:
            m = types_s == c
            ax.scatter(coords_s[m, 0], coords_s[m, 1], s=1, c=[colors_map[c]], label=c, rasterized=True)
        ax.set_title(f"Spatial cell types ({method})")
        ax.legend(fontsize=6, markerscale=4, loc="upper right")
        ax.set_aspect("equal")
        ax.axis("off")
        plt.tight_layout()
        png2 = os.path.join(args.out_dir, "spatial_celltypes.png")
        plt.savefig(png2, dpi=150, bbox_inches="tight")
        plt.close(fig)

    print("ALL_DONE")


if __name__ == "__main__":
    main()
