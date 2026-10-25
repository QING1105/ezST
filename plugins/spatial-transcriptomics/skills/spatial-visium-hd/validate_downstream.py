#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""validate_downstream.py — Visium HD 细胞级结果下游验证。

作用
----
Visium HD 分支（spatial-visium-hd skill）在 S4 输出细胞级 h5ad 后，
需要验证该结果能正常进入共享下游流程（normalize -> cluster），
确认细胞矩阵质量可用后再交给用户审查。

本脚本完成两件事：
  1. 从 S4 输出的细胞级 h5ad 取一个子集（默认 1 万细胞，速度快）；
  2. 调用 ezST 共享下游函数 normalize_visium + cluster_spatial_data，
     生成 normalized/clustered h5ad 和 UMAP/HVG 图，验证流程通。

用法
----
  python validate_downstream.py \
      --input <S4 细胞级 h5ad> \
      --out-dir <输出目录> \
      --n-cells 10000 \
      --n-top-genes 2000 \
      --n-pcs 30 \
      --resolution 1.0

输出
----
  <out-dir>/sub10k.h5ad                # 子集 (原始计数)
  <out-dir>/sub10k_normalized.h5ad     # normalize + log1p 后
  <out-dir>/sub10k_hvg.png             # 高变基因图
  <out-dir>/sub10k_clustered.h5ad      # PCA + Leiden 聚类后
  <out-dir>/sub10k_umap.png            # UMAP 图
"""

import argparse
import os
import sys
import time

try:
    from .visium import normalize_visium, cluster_spatial_data
except ImportError:
    from visium import normalize_visium, cluster_spatial_data


def parse_args():
    """解析命令行参数。"""
    ap = argparse.ArgumentParser(description="Visium HD 下游验证")
    ap.add_argument("--input", required=True,
                    help="S4 输出的细胞级 h5ad（04_bin2cell.h5ad / 04_proseg_cells.h5ad）")
    ap.add_argument("--out-dir", default="downstream_test",
                    help="输出目录（默认 downstream_test）")
    ap.add_argument("--n-cells", type=int, default=10000,
                    help="验证用细胞子集大小（默认 10000）")
    ap.add_argument("--n-top-genes", type=int, default=2000,
                    help="高变基因数（默认 2000）")
    ap.add_argument("--n-pcs", type=int, default=30,
                    help="PCA 组件数（默认 30）")
    ap.add_argument("--resolution", type=float, default=None,
                    help="Leiden 分辨率（默认 None=1.0, LLM 可覆盖）")
    return ap.parse_args()


def main():
    """下游验证主流程。"""
    args = parse_args()
    t0 = time.time()

    # --- 1. 读取 S4 细胞级结果 ---
    print(f"[1/5] 读取细胞级结果: {args.input}", flush=True)
    import anndata as ad
    import numpy as np

    adata = ad.read_h5ad(args.input)
    print(f"  完整矩阵: {adata.shape[0]} cells x {adata.shape[1]} genes", flush=True)

    # --- 2. 取子集（保证验证速度快；全量跑留给正式下游） ---
    print(f"[2/5] 取 {args.n_cells} 细胞子集", flush=True)
    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(42)  # 固定随机种子，结果可复现
    idx = rng.choice(adata.n_obs, min(args.n_cells, adata.n_obs), replace=False)
    sub = adata[idx].copy()
    sub_path = os.path.join(args.out_dir, f"sub{args.n_cells}.h5ad")
    sub.write(sub_path)
    print(f"  子集: {sub.shape} -> {sub_path} ({time.time()-t0:.0f}s)", flush=True)

    # --- 3. normalize（target_sum=10000, log1p, 选 HVG） ---
    print(f"[3/5] normalize (top {args.n_top_genes} HVG)", flush=True)
    norm_path = os.path.join(args.out_dir, f"sub{args.n_cells}_normalized.h5ad")
    norm_plot = os.path.join(args.out_dir, f"sub{args.n_cells}_hvg.png")
    msg = normalize_visium(
        adata_path=sub_path,
        output_path=norm_path,
        plot_path=norm_plot,
        n_top_genes=args.n_top_genes,
    )
    print(f"  normalize OK: {msg[:120]} ({time.time()-t0:.0f}s)", flush=True)

    # --- 4. cluster（PCA + UMAP + Leiden） ---
    print(f"[4/5] cluster (PCA {args.n_pcs}, res={args.resolution})", flush=True)
    clu_path = os.path.join(args.out_dir, f"sub{args.n_cells}_clustered.h5ad")
    clu_umap = os.path.join(args.out_dir, f"sub{args.n_cells}_umap.png")
    msg2 = cluster_spatial_data(
        adata_path=norm_path,
        output_path=clu_path,
        plot_path=clu_umap,
        n_pcs=args.n_pcs,
        resolution=args.resolution,
    )
    print(f"  cluster OK: {msg2[:120]} ({time.time()-t0:.0f}s)", flush=True)

    # --- 5. 校验输出完整性 ---
    print("[5/5] 校验输出", flush=True)
    outputs = [sub_path, norm_path, norm_plot, clu_path, clu_umap]
    ok = True
    for f in outputs:
        exists = os.path.exists(f)
        print(f"  {'OK ' if exists else 'MISSING'} {f}", flush=True)
        ok = ok and exists

    if not ok:
        print("DOWNSTREAM VALIDATION FAILED: 输出缺失", flush=True)
        sys.exit(1)

    print(f"DOWNSTREAM VALIDATION DONE ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
