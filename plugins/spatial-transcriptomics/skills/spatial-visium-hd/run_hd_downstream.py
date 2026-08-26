#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""run_hd_downstream.py — Visium HD 分割后完整下游

S1-S5 (分割+QC) 之后:
  1. 注释: label transfer (Haber 2017 参考, scVI + LogisticRegression)
  2. domains: squidpy 空间邻域图 + Leiden
  3. neighborhood: squidpy neighborhood enrichment (细胞类型共定位)
  4. CellChat: 细胞类型通讯 (R, spatial 模式, 自动 downsample)

依据: 分割后为"空间单细胞", 不走 deconvolve (见概念文档)
"""
import argparse
import os
import subprocess
import sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True, help="Visium HD 细胞级 h5ad")
    ap.add_argument("--ref", default="", help="参考 h5ad (label transfer 用), 必填除非 --skip-annotate")
    ap.add_argument("--cell-type-key", default="cell_type")
    ap.add_argument("--out-dir", default="./hd_downstream")
    ap.add_argument("--n-neighs", type=int, default=None,
                    help="空间邻域 K (None=数据驱动, LLM 可覆盖)")
    ap.add_argument("--resolution", type=float, default=None, help="domain resolution (LLM decision)")
    ap.add_argument("--skip-annotate", action="store_true")
    ap.add_argument("--skip-domains", action="store_true")
    ap.add_argument("--skip-neighbor", action="store_true")
    ap.add_argument("--skip-cellchat", action="store_true")
    ap.add_argument("--species", default=None,
                    help="物种 human/mouse (None=自动检测, 检测失败需 LLM 传参)")
    ap.add_argument("--rscript", default="Rscript", help="Rscript 可执行文件路径 (CellChat 用, 默认走 PATH)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    import numpy as np
    import scanpy as sc
    import anndata as ad

    # ============ 0. 加载 ============
    print("=== Load cells ===")
    adata = ad.read_h5ad(args.cells)
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes")
    if "spatial" not in adata.obsm:
        print("ERROR: no obsm['spatial']"); sys.exit(1)

    # ============ 1. 注释 (label transfer, 双方法 auto) ============
    if not args.skip_annotate and "cell_type" not in adata.obs.columns:
        if not args.ref:
            print("ERROR: --ref 必填 (注释需要参考数据), 或加 --skip-annotate 用已有注释")
            sys.exit(1)
        print("=== 1. Annotation (dual-method auto: scVI if counts ref, ingest if log1p ref) ===")
        # 复用 annotate_cells.py 双方法逻辑 (scVI / PCA+kNN ingest, --method auto)
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "annotate_cells.py")
        subprocess.run(
            [sys.executable, script,
             "--cells", args.cells,
             "--ref", args.ref,
             "--cell-type-key", args.cell_type_key,
             "--method", "auto",
             "--out-dir", args.out_dir],
            check=True,
        )
        ann = ad.read_h5ad(os.path.join(args.out_dir, "annotated_cells.h5ad"))
        adata.obs["cell_type"] = ann.obs["cell_type"]
        adata.obs["cell_type_prob"] = ann.obs["cell_type_prob"]
        vc = adata.obs["cell_type"].value_counts()
        print("  top types:")
        print(vc.head(12).to_string())
        adata.write(os.path.join(args.out_dir, "01_annotated.h5ad"), compression="lzf")
        print("  saved 01_annotated.h5ad")

    # ============ 2. Domains (统一调用 skill 的 identify_spatial_domains) ============
    if not args.skip_domains:
        print("=== 2. Spatial domains ===")
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from visium import identify_spatial_domains
        # 传路径 (函数内部读 h5ad -> 分析 -> 写回 + 画图)
        msg = identify_spatial_domains(
            args.cells,
            output_path=os.path.join(args.out_dir, "02_domains.h5ad"),
            plot_path=os.path.join(args.out_dir, "02_domains.png"),
            resolution=args.resolution,
        )
        print(msg)
        n_dom = adata.obs["spatial_domain"].nunique()
        print(f"  domains: {n_dom}")

    # ============ 3. Neighborhood enrichment ============
    if not args.skip_neighbor and "cell_type" in adata.obs.columns:
        print("=== 3. Neighborhood enrichment ===")
        import squidpy as sq
        if "spatial_neighbors" not in adata.uns:
            sq.gr.spatial_neighbors(adata, n_neighs=args.n_neighs, coord_type="generic")
        sq.gr.nhood_enrichment(adata, cluster_key="cell_type")
        # squidpy >=1.6 存到 uns['cell_type_nhood_enrichment']['zscore'] (可能是 ndarray)
        import numpy as np
        import pandas as pd
        if "cell_type_nhood_enrichment" in adata.uns:
            z = adata.uns["cell_type_nhood_enrichment"]["zscore"]
        else:
            z = adata.uns["neighborhood_enrichment"]["zscore"]
        cats = list(adata.obs["cell_type"].astype("category").cat.categories)
        if isinstance(z, np.ndarray):
            df = pd.DataFrame(z, index=cats, columns=cats)
        else:
            df = z
        df.to_csv(os.path.join(args.out_dir, "03_neighborhood_enrichment.csv"))
        import matplotlib.pyplot as plt
        import seaborn as sns
        fig, ax = plt.subplots(figsize=(max(8, df.shape[0]*0.5), max(8, df.shape[0]*0.5)))
        sns.heatmap(df, cmap="RdBu_r", center=0, ax=ax, square=True)
        ax.set_title("Neighborhood enrichment (z-score)")
        plt.tight_layout()
        plt.savefig(os.path.join(args.out_dir, "03_neighborhood.png"), dpi=150)
        plt.close(fig)
        print("  saved 03_neighborhood_enrichment.csv + 03_neighborhood.png")

    # ============ 4. CellChat ============
    if not args.skip_cellchat and "cell_type" in adata.obs.columns:
        print("=== 4. CellChat (spatial, downsample) ===")
        # 每类型最多 1000 细胞 (官方建议 >50k 时 downsample)
        downsample = adata.obs.groupby("cell_type", observed=True).apply(
            lambda x: x.sample(min(len(x), 1000), random_state=0)).index.get_level_values(-1)
        cc = adata[downsample].copy()
        print(f"  CellChat input: {cc.n_obs} cells ({len(cc.obs['cell_type'].unique())} types)")
        cc.write(os.path.join(args.out_dir, "04_cellchat_input.h5ad"), compression="lzf")
        # 调用经典 visium.py 的 CellChat 桥
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from visium import infer_spatial_cell_communication
        msg = infer_spatial_cell_communication(
            os.path.join(args.out_dir, "04_cellchat_input.h5ad"),
            cell_type_key="cell_type",
            output_dir=os.path.join(args.out_dir, "cellchat_results"),
            species=args.species,  # None -> 自动检测 (LLM-in-the-loop)
            r_script_path=args.rscript,
        )
        print(msg)

    print("=== DONE ===")

if __name__ == "__main__":
    main()
