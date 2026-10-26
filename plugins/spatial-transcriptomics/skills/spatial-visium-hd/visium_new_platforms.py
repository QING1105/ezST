"""
Spatial-platform loading, segmentation, aggregation, and visualization helpers.
"""

import os

try:
    from ._geometry import (load_nuclei_polygons, save_nuclei_polygons,
                            filter_nuclei_instances, nuclei_xy_polygons, read_proseg_counts)
except ImportError:  # direct execution of the bundled skill scripts
    from _geometry import (load_nuclei_polygons, save_nuclei_polygons,
                           filter_nuclei_instances, nuclei_xy_polygons, read_proseg_counts)


def _probe_gpu() -> bool:
    """探测 TensorFlow GPU 是否可用且能编译卷积。

    必须在 TF 首次初始化之前调用。服务器上若 cuDNN 与驱动不兼容,
    GPU 编译卷积会抛 Autotuner error, 此时应强制 CPU。
    返回 True = 可以用 GPU; False = 需要 CPU。

    注意: 必须用 tf.function(jit_compile=True) 强制 XLA 编译路径,
    Keras model.predict 走的正是这条路径。eager 小卷积不触发 XLA,
    会漏掉 cuDNN 编译失败。
    """
    # cuDNN frontend 在 NCHW 布局下对 1x1 卷积找不到配置 (Autotuner error),
    # 强制 NHWC 布局绕过 (TF 2.21 + cuDNN 9.5 实测)。提前设好, 避免 TF
    # 初始化后 XLA_FLAGS 改动不生效。
    import os as _os
    _xla = _os.environ.get("XLA_FLAGS", "")
    if "force_conv_nhwc" not in _xla:
        _xla = (_xla + " --xla_gpu_force_conv_nhwc=true").strip()
        _os.environ["XLA_FLAGS"] = _xla
    try:
        import tensorflow as _tf
        if not _tf.config.list_physical_devices("GPU"):
            return False
        # 模拟 StarDist 的 3x3 conv + bias + ReLU (256 通道), 强制 XLA 编译
        _conv = _tf.keras.layers.Conv2D(
            256, 3, padding="same", activation="relu"
        )

        @_tf.function(jit_compile=True)
        def _probe_conv(x):
            return _conv(x)

        _x = _tf.random.normal((1, 16, 16, 256))
        _ = _probe_conv(_x)
        _tf.keras.backend.clear_session()
        return True
    except Exception:
        try:
            _tf.keras.backend.clear_session()
        except Exception:
            pass
        return False


def _plot_setup(output_dir: str) -> None:
    """Ensure matplotlib uses the Agg backend and the output dir exists."""
    import matplotlib

    matplotlib.use("Agg")
    os.makedirs(output_dir, exist_ok=True)


def _plot_cells_overlay(adata, image_path, plot_path, title=None, downsample=1,
                        max_points=None, geometry_path=None, geometry_kind="auto",
                        edge_color="tab:red", linewidth=0.12, dpi=150,
                        figsize_scale=1.0, allow_legacy_pickle=False):
    """HE background + cell overlay plot (HD, full-resolution).

    Supports two visual modes:
      1. Geometry mode (default when geometry_path given): draws the true cell
         boundaries as polygons (Proseg WKB parquet or StarDist polys.npz).
      2. Scatter fallback: centroids from adata.obsm['spatial'].

    Args:
        adata: AnnData with obsm['spatial'] (fullres pixel coords).
        image_path: Path to H&E image (btf/tif/png/jpg).
        plot_path: Where to save the PNG.
        title: Optional title (default '<n> cells').
        downsample: Downsample factor for the HE background (1 = full res).
        max_points: Cap on scatter points (None = all).
        geometry_path: Path to geometry source:
            - Proseg: proseg_out.zarr/shapes/cell_boundaries/shapes.parquet
            - StarDist: polys.npz (list of shapely polygons / dict with 'coord')
        geometry_kind: 'proseg' | 'stardist' | 'auto' (infer from extension).
        edge_color: Polygon edge color.
        linewidth: Polygon edge line width.
        dpi: Output PNG dpi.
        figsize_scale: Multiplier for figure size (default 1.0 -> 16x16 at full res).
    """
    _plot_setup(os.path.dirname(plot_path) or ".")
    import numpy as np
    import matplotlib.pyplot as plt

    # HE 背景 (fullres)
    img = None
    try:
        ext = os.path.splitext(image_path)[1].lower()
        if ext in (".tif", ".tiff", ".btf"):
            from tifffile import imread
            img = imread(image_path)
            if img.ndim == 3 and img.shape[-1] == 4:
                img = img[..., :3]
        else:
            import imageio.v2 as imageio
            img = imageio.imread(image_path)
            if img.ndim == 3 and img.shape[-1] == 4:
                img = img[..., :3]
        if img.ndim == 2:
            img = np.stack([img] * 3, axis=-1)
        if img.shape[0] > 0 and img.shape[1] > 0 and downsample > 1:
            img = img[::downsample, ::downsample]
    except Exception:
        img = None  # 读图失败就纯几何/散点

    polys = None
    if geometry_path and os.path.exists(geometry_path):
        try:
            kind = geometry_kind
            if kind == "auto":
                kind = "proseg" if geometry_path.endswith(".parquet") else "stardist"
            if kind == "proseg":
                import pandas as pd
                from shapely import wkb
                df = pd.read_parquet(geometry_path)
                scale = getattr(adata, "uns", {}).get("cellpose_scale", None)
                if scale is None:
                    scale = 0.2738
                polys = []
                for b in df["geometry"].values:
                    g = wkb.loads(b)
                    if g.geom_type == "Polygon":
                        polys.append(np.array(g.exterior.coords) / scale)
                    elif g.geom_type == "MultiPolygon":
                        parts = sorted(g.geoms, key=lambda p: p.area, reverse=True)
                        polys.append(np.array(parts[0].exterior.coords) / scale)
            else:
                polys = nuclei_xy_polygons(geometry_path, allow_legacy_pickle=allow_legacy_pickle,
                                           object_ids=adata.obs.get("object_id"))
        except Exception as exc:
            raise ValueError(f"Cannot load geometry {geometry_path}: {exc}") from exc

    fig_w = 20 * figsize_scale
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_w), dpi=dpi)
    if img is not None:
        ax.imshow(img, cmap="gray")

    if polys is not None and len(polys) > 0:
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Polygon as MplPolygon
        ds = max(downsample, 1)
        patches = [MplPolygon(c / ds, closed=True) for c in polys]
        pc = PatchCollection(patches, facecolor="none", edgecolor=edge_color,
                             linewidths=linewidth, alpha=0.8)
        ax.add_collection(pc)
    else:
        coords = adata.obsm["spatial"]
        if max_points and len(coords) > max_points:
            rng = np.random.default_rng(0)
            idx = rng.choice(len(coords), max_points, replace=False)
            coords = coords[idx]
        ax.scatter(coords[:, 0] / max(downsample, 1),
                   coords[:, 1] / max(downsample, 1),
                   s=0.3, alpha=0.6, c=edge_color, rasterized=True)

    if title is None:
        title = f"{adata.n_obs} cells"
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_cells_crop(adata, image_path, plot_path, roi=None, roi_size=1000,
                    geometry_path=None, geometry_kind="auto",
                    edge_color="tab:blue", linewidth=0.8, dpi=150, title=None,
                    allow_legacy_pickle=False):
    """HE background + cell overlay for a user-selected ROI crop (HD).

    Lets the user choose a region of interest (in fullres pixel coords) and
    plots only that crop with the true cell boundaries (or centroids).

    Args:
        adata: AnnData with obsm['spatial'] (fullres pixel coords).
        image_path: Path to H&E image (btf/tif/png/jpg).
        plot_path: Where to save the PNG.
        roi: Region of interest in fullres pixel coords:
            - tuple of 4: (x0, y0, x1, y1) explicit box
            - tuple of 2: (cx, cy) center -> expand by roi_size
            - None: auto-pick the densest region from cell centroids
        roi_size: Square crop size when roi is (cx, cy) or None (default 1000).
        geometry_path: Path to geometry source (Proseg parquet / StarDist NPZ).
        geometry_kind: 'proseg' | 'stardist' | 'auto'.
        edge_color: Polygon edge color.
        linewidth: Polygon edge line width.
        dpi: Output PNG dpi.
        title: Optional title (default includes ROI box and cell count).
    """
    _plot_setup(os.path.dirname(plot_path) or ".")
    import numpy as np
    import matplotlib.pyplot as plt

    # HE 背景
    img = None
    try:
        ext = os.path.splitext(image_path)[1].lower()
        if ext in (".tif", ".tiff", ".btf"):
            from tifffile import imread
            img = imread(image_path)
            if img.ndim == 3 and img.shape[-1] == 4:
                img = img[..., :3]
        else:
            import imageio.v2 as imageio
            img = imageio.imread(image_path)
            if img.ndim == 3 and img.shape[-1] == 4:
                img = img[..., :3]
        if img.ndim == 2:
            img = np.stack([img] * 3, axis=-1)
    except Exception:
        img = None

    coords = adata.obsm["spatial"]
    h_img, w_img = (img.shape[0], img.shape[1]) if img is not None else (0, 0)

    # --- 确定 ROI ---
    if roi is None:
        # 自动选细胞最密集区
        if len(coords) > 1000:
            rng = np.random.default_rng(0)
            samp = coords[rng.choice(len(coords), min(len(coords), 20000), replace=False)]
            hh, xedges, yedges = np.histogram2d(samp[:, 0], samp[:, 1], bins=20)
            idx = np.unravel_index(np.argmax(hh), hh.shape)
            cx = (xedges[idx[0]] + xedges[idx[0] + 1]) / 2
            cy = (yedges[idx[1]] + yedges[idx[1] + 1]) / 2
        else:
            cx, cy = coords.mean(axis=0)
        s = roi_size // 2
        x0, y0 = int(cx - s), int(cy - s)
        x1, y1 = int(cx + s), int(cy + s)
    elif len(roi) == 4:
        x0, y0, x1, y1 = [int(v) for v in roi]
    else:  # (cx, cy)
        cx, cy = int(roi[0]), int(roi[1])
        s = roi_size // 2
        x0, y0, x1, y1 = cx - s, cy - s, cx + s, cy + s

    if img is not None:
        x0, y0 = max(0, x0), max(0, y0)
        x1, y1 = min(w_img, x1), min(h_img, y1)
    if x1 <= x0 or y1 <= y0:
        raise ValueError(f"bad ROI: x[{x0}:{x1}] y[{y0}:{y1}]")

    # --- 几何 (ROI 内) ---
    polys = []
    if geometry_path and os.path.exists(geometry_path):
        try:
            kind = geometry_kind
            if kind == "auto":
                kind = "proseg" if geometry_path.endswith(".parquet") else "stardist"
            if kind == "proseg":
                import pandas as pd
                from shapely import wkb
                df = pd.read_parquet(geometry_path)
                scale = getattr(adata, "uns", {}).get("cellpose_scale", None)
                if scale is None:
                    scale = 0.2738
                for b in df["geometry"].values:
                    g = wkb.loads(b)
                    if g.geom_type == "Polygon":
                        c = np.array(g.exterior.coords) / scale
                    elif g.geom_type == "MultiPolygon":
                        parts = sorted(g.geoms, key=lambda p: p.area, reverse=True)
                        c = np.array(parts[0].exterior.coords) / scale
                    else:
                        continue
                    cx_, cy_ = c[:, 0].mean(), c[:, 1].mean()
                    if x0 <= cx_ < x1 and y0 <= cy_ < y1:
                        polys.append(c)
            else:
                for c in nuclei_xy_polygons(geometry_path, allow_legacy_pickle=allow_legacy_pickle,
                                           object_ids=adata.obs.get("object_id")):
                    # Draw polygons intersecting the ROI, including cells centered just outside it.
                    if c[:, 0].max() >= x0 and c[:, 0].min() < x1 and c[:, 1].max() >= y0 and c[:, 1].min() < y1:
                        polys.append(c)
        except Exception as exc:
            raise ValueError(f"Cannot load geometry {geometry_path}: {exc}") from exc

    # --- 绘制 ---
    fig, ax = plt.subplots(1, 1, figsize=(12, 12), dpi=dpi)
    if img is not None:
        ax.imshow(img[y0:y1, x0:x1])
    if polys:
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Polygon as MplPolygon
        patches = [MplPolygon(c - [x0, y0], closed=True) for c in polys]
        pc = PatchCollection(patches, facecolor="none", edgecolor=edge_color,
                             linewidths=linewidth, alpha=0.9)
        ax.add_collection(pc)
        n_cells = len(polys)
    else:
        mask = (coords[:, 0] >= x0) & (coords[:, 0] < x1) & \
               (coords[:, 1] >= y0) & (coords[:, 1] < y1)
        cc = coords[mask]
        ax.scatter(cc[:, 0] - x0, cc[:, 1] - y0, s=2, alpha=0.7,
                   c=edge_color, linewidths=0, rasterized=True)
        n_cells = len(cc)

    if title is None:
        title = f"{n_cells} cells (crop x{x0}:{x1}, y{y0}:{y1})"
    ax.set_xlim(-0.5, x1 - x0 - 0.5)
    ax.set_ylim(y1 - y0 - 0.5, -0.5)
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _plot_nuclei_crop(labels, image_path, plot_path, roi_size=512,
                      title="Nuclei segmentation"):
    """Nuclei segmentation check plot: HE crop + nuclei boundaries overlay."""
    _plot_setup(os.path.dirname(plot_path) or ".")
    import numpy as np
    import matplotlib.pyplot as plt
    from skimage.segmentation import find_boundaries

    img = None
    try:
        ext = os.path.splitext(image_path)[1].lower()
        if ext in (".tif", ".tiff", ".btf"):
            from tifffile import imread
            img = imread(image_path)
            if img.ndim == 3 and img.shape[-1] == 4:
                img = img[..., :3]
        else:
            import imageio.v2 as imageio
            img = imageio.imread(image_path)
            if img.ndim == 3 and img.shape[-1] == 4:
                img = img[..., :3]
        if img.ndim == 2:
            img = np.stack([img] * 3, axis=-1)
    except Exception:
        img = None

    # 中心 ROI
    h, w = labels.shape
    cy, cx = h // 2, w // 2
    y0, y1 = max(0, cy - roi_size // 2), min(h, cy + roi_size // 2)
    x0, x1 = max(0, cx - roi_size // 2), min(w, cx + roi_size // 2)
    lab_crop = labels[y0:y1, x0:x1]
    bounds = find_boundaries(lab_crop, mode="outer")

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    if img is not None:
        ax.imshow(img[y0:y1, x0:x1])
    ax.contour(bounds, levels=[0.5], colors="lime", linewidths=0.6)
    ax.set_title(f"{title} (crop {x0}:{x1}, {y0}:{y1})")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def plot_cell_qc(
    adata,
    output_dir: str,
    sample: str = "sample",
    method: str = "stardist",
) -> str:
    """Cell-level QC visualization for Visium HD segmentation output.

    Follows community standard QC practice (SIB-Swiss spatial training,
    bcbio spatial-reports, SpatialQC):
      - Panel A: violin of n_genes_by_counts (genes per cell)
      - Panel B: violin of total_counts (UMI per cell)
      - Panel C: violin of pct_counts_mt (mitochondrial %, if mt genes present)
      - Panel D: histogram of novelty score (n_genes/n_UMI, good > 0.8)
      - Panel E: histogram of bin_count (bins assigned per cell)
      - Panel F: spatial scatter of low-quality cells (flagged), vs all cells
    Also writes ``qc_metrics.json`` with summary numbers and flags
    (low-count threshold, mito threshold, doublet-like cells).

    Args:
        adata: Cell-level AnnData (from bin_to_cell / run_proseg).
        output_dir: Where to save ``05_qc.png`` and ``qc_metrics.json``.
        sample: Sample identifier for titles.
        method: 'stardist' or 'proseg' for titles.

    Returns:
        str: Log message.
    """
    _plot_setup(output_dir)
    import numpy as np
    import json
    import matplotlib.pyplot as plt

    # --- 计算 QC 指标 ---
    import scanpy as sc
    adata.var_names_make_unique()
    # n_genes / total_counts / mito% (若尚未存在)
    if "n_genes_by_counts" not in adata.obs.columns:
        mt_genes = [g for g in adata.var_names if g.startswith("mt-")]
        if mt_genes:
            adata.var["mt"] = adata.var_names.isin(mt_genes)
            sc.pp.calculate_qc_metrics(adata, inplace=True, qc_vars=["mt"])
        else:
            sc.pp.calculate_qc_metrics(adata, inplace=True)

    n_cells = adata.n_obs
    n_genes = adata.obs["n_genes_by_counts"].astype(float)
    n_umis = adata.obs["total_counts"].astype(float)
    novelty = (n_genes / n_umis).replace([np.inf, np.nan], 0.0)
    has_mt = "pct_counts_mt" in adata.obs.columns
    mito = adata.obs["pct_counts_mt"].astype(float) if has_mt else None

    # --- 阈值 (bcbio/SpatialQC 惯例) ---
    # novelty > 0.8 为好; mito < 20% 为好; bins >= 5 为细胞 (10x 默认扩展 8um ~ 4 bins)
    low_novelty = novelty < 0.8
    high_mito = (mito > 20.0) if has_mt else np.zeros(n_cells, dtype=bool)
    low_bins = adata.obs["bin_count"] < 5 if "bin_count" in adata.obs.columns else np.zeros(n_cells, dtype=bool)
    low_genes = n_genes < 200  # scRNA 常见过滤阈值
    flagged = low_novelty | high_mito | low_bins | low_genes
    n_flagged = int(flagged.sum())

    med_genes = float(np.median(n_genes))
    med_counts = float(np.median(n_umis))
    med_novelty = float(np.median(novelty))
    med_mito = float(np.median(mito)) if has_mt else None
    bin_col = "bin_count" if "bin_count" in adata.obs.columns else None
    med_bins = float(np.median(adata.obs[bin_col])) if bin_col else None

    metrics = {
        "sample": sample,
        "method": method,
        "n_cells": int(n_cells),
        "median_genes_per_cell": med_genes,
        "median_counts_per_cell": med_counts,
        "median_novelty": med_novelty,
        "median_mito_pct": med_mito,
        "median_bins_per_cell": med_bins,
        "n_flagged_low_quality": n_flagged,
        "pct_flagged": round(100.0 * n_flagged / max(n_cells, 1), 2),
        "thresholds": {
            "novelty_min": 0.8,
            "mito_max_pct": 20.0,
            "bins_min": 5,
            "genes_min": 200,
        },
        "note": "10x Space Ranger default nucleus expansion = 8 um (~4 bins); "
                "median bins/cell < 5 suggests missing nucleus expansion in bin_to_cell",
    }
    with open(os.path.join(output_dir, "qc_metrics.json"), "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2, ensure_ascii=False)

    # --- 2x3 QC 面板 ---
    fig, axes = plt.subplots(2, 3, figsize=(18, 11), dpi=150)

    # A: n_genes violin
    ax = axes[0, 0]
    ax.violinplot(n_genes.values, showmedians=True)
    ax.axhline(med_genes, color="red", ls="--", lw=1, label=f"median={med_genes:.0f}")
    ax.axhline(200, color="orange", ls=":", lw=1, label="min 200")
    ax.set_title(f"Genes per cell (median {med_genes:.0f})")
    ax.set_ylabel("n_genes")
    ax.set_xticks([])
    ax.legend(fontsize=8)

    # B: total_counts violin
    ax = axes[0, 1]
    ax.violinplot(n_umis.values, showmedians=True)
    ax.axhline(med_counts, color="red", ls="--", lw=1, label=f"median={med_counts:.0f}")
    ax.set_title(f"UMI per cell (median {med_counts:.0f})")
    ax.set_ylabel("total_counts")
    ax.set_xticks([])
    ax.legend(fontsize=8)

    # C: mito% violin (if available)
    ax = axes[0, 2]
    if has_mt:
        ax.violinplot(mito.values, showmedians=True)
        ax.axhline(med_mito, color="red", ls="--", lw=1, label=f"median={med_mito:.1f}%")
        ax.axhline(20, color="orange", ls=":", lw=1, label="max 20%")
        ax.set_title(f"Mito % (median {med_mito:.1f}%)")
        ax.set_ylabel("pct_counts_mt")
    else:
        ax.text(0.5, 0.5, "no mt- genes", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Mito % (N/A)")
    ax.set_xticks([])
    ax.legend(fontsize=8)

    # D: novelty score histogram
    ax = axes[1, 0]
    ax.hist(novelty.values, bins=60, color="#4C72B0", alpha=0.8)
    ax.axvline(med_novelty, color="red", ls="--", lw=1, label=f"median={med_novelty:.2f}")
    ax.axvline(0.8, color="orange", ls=":", lw=1, label="good > 0.8")
    ax.set_title(f"Novelty score (n_genes/n_UMI, median {med_novelty:.2f})")
    ax.set_xlabel("novelty")
    ax.set_ylabel("cells")
    ax.legend(fontsize=8)

    # E: bin_count histogram
    ax = axes[1, 1]
    if bin_col:
        bc = adata.obs[bin_col].astype(float)
        ax.hist(bc.values, bins=60, color="#C44E52", alpha=0.8)
        if med_bins is not None:
            ax.axvline(med_bins, color="red", ls="--", lw=1, label=f"median={med_bins:.0f}")
        ax.axvline(5, color="orange", ls=":", lw=1, label="min 5 bins")
        ax.set_title(f"Bins per cell (median {med_bins:.0f})")
        ax.set_xlabel("bin_count")
        ax.set_ylabel("cells")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "no bin_count", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Bins per cell (N/A)")
        ax.axis("off")

    # F: spatial low-quality flag
    ax = axes[1, 2]
    coords = adata.obsm["spatial"]
    if len(coords) > 80000:
        rng = np.random.default_rng(0)
        idx = rng.choice(len(coords), 80000, replace=False)
    else:
        idx = np.arange(len(coords))
    ax.scatter(coords[idx, 0], coords[idx, 1], s=1, alpha=0.5,
               c="lightgray", rasterized=True, label=f"all ({len(idx):,})")
    flag_idx = idx[flagged[idx]]
    if len(flag_idx):
        ax.scatter(coords[flag_idx, 0], coords[flag_idx, 1], s=3, alpha=0.8,
                   c="red", rasterized=True, label=f"flagged ({n_flagged:,})")
    ax.set_title(f"Spatial low-quality cells ({n_flagged:,}, {100.0*n_flagged/max(n_cells,1):.1f}%)")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(fontsize=8, loc="upper right")

    fig.suptitle(f"{sample} [{method}] QC — {n_cells:,} cells", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    qc_png = os.path.join(output_dir, "05_qc.png")
    plt.savefig(qc_png, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return (f"QC panel saved -> {qc_png} "
            f"(n_cells={n_cells}, median_genes={med_genes:.0f}, "
            f"median_counts={med_counts:.0f}, median_novelty={med_novelty:.2f}, "
            f"flagged={n_flagged:,} [{100.0*n_flagged/max(n_cells,1):.1f}%])")


def load_visium_hd(
    h5_path: str,
    tissue_positions_path: str = None,
    sample_id: str = None,
    output_path: str = "./visium_hd_loaded.h5ad",
    plot_path: str = None,
) -> str:
    """Load Visium HD data from the 2 um binned matrix.

    Visium HD Space Ranger outputs a 2 um bin x gene matrix. The 2 um bin is a
    transcript-capture unit, NOT a biological unit. This function loads the
    binned matrix and attaches bin coordinates; cell reconstruction
    (segmentation + bin-to-cell) is done by :func:`segment_cells` +
    :func:`bin_to_cell`.

    Supported inputs:
      1. h5ad with obsm['spatial'] already attached.
      2. 10x h5 (filtered_feature_bc_matrix.h5) + tissue_positions.parquet.

    Args:
        h5_path: Path to .h5ad or filtered_feature_bc_matrix.h5.
        tissue_positions_path: Optional path to tissue_positions.parquet
            (Visium HD bin coordinates; barcode, x, y columns).
        sample_id: Optional sample identifier.
        output_path: Where to save the loaded AnnData.
        plot_path: Optional layout plot path.

    Returns:
        str: Step-by-step log.
    """
    import scanpy as sc
    import anndata as ad

    steps = []
    if not os.path.exists(h5_path):
        return f"Error: h5_path not found: {h5_path}"

    if h5_path.endswith(".h5ad"):
        adata = ad.read_h5ad(h5_path)
        steps.append(f"Loaded h5ad: {adata.shape[0]} bins x {adata.shape[1]} genes")
    elif h5_path.endswith(".h5"):
        adata = sc.read_10x_h5(h5_path)
        adata.var_names_make_unique()
        steps.append(f"Loaded 10x h5: {adata.shape[0]} bins x {adata.shape[1]} genes")
    else:
        return f"Error: unsupported Visium HD input format: {h5_path}"

    # Attach coordinates from parquet if provided
    if tissue_positions_path and os.path.exists(tissue_positions_path):
        import pandas as pd
        pos = pd.read_parquet(tissue_positions_path)
        # 语义匹配列名: barcode + x/y (pxl_col/pxl_row 或 x/y 或 col/row)
        col_map = {}
        for i, c in enumerate(pos.columns):
            cl = str(c).strip().lower()
            if cl in ("barcode", "bc", "barcode_id", "cell_id", "id"):
                col_map["barcode"] = c
            elif cl in ("pxl_col_in_fullres", "pxl_col", "x_pixel", "x"):
                col_map["x"] = c
            elif cl in ("pxl_row_in_fullres", "pxl_row", "y_pixel", "y"):
                col_map["y"] = c
        if "x" in col_map and "y" in col_map:
            if "barcode" in col_map:
                pos.index = pos[col_map["barcode"]].astype(str).values
                coords = pos.loc[adata.obs_names, [col_map["x"], col_map["y"]]]
                adata.obsm["spatial"] = coords.astype(float).values
            else:
                if len(pos) == adata.n_obs:
                    adata.obsm["spatial"] = pos[[col_map["x"], col_map["y"]]].astype(float).values
                else:
                    steps.append(f"Warning: parquet has {len(pos)} rows but adata has {adata.n_obs} bins; coords NOT attached")
            steps.append(f"Attached coordinates from {os.path.basename(tissue_positions_path)}")
        else:
            steps.append(f"Warning: no recognizable x/y columns in parquet: {list(pos.columns)}")

    if sample_id:
        adata.obs["sample_id"] = sample_id

    adata.uns["spatial_meta"] = {
        "source": h5_path,
        "platform": "visium_hd",
        "n_bins": adata.n_obs,
        "n_genes": adata.n_vars,
        "has_coords": "spatial" in adata.obsm,
    }

    adata.write(output_path, compression="lzf")
    steps.append(f"Saved loaded Visium HD AnnData to {output_path}")
    return "\n".join(steps)


def _load_stardist_pretrained(model_class, model_type):
    """Use Keras' checksum-verified cache when Windows cannot create symlinks."""
    try:
        return model_class.from_pretrained(model_type)
    except OSError as error:
        if getattr(error, "winerror", None) != 1314:
            raise
        from pathlib import Path
        from csbdeep.models.pretrained import get_model_details, get_file
        key, _, details = get_model_details(model_class, model_type)
        target = str(Path("models") / model_class.__name__ / key)
        extracted = Path(get_file(fname=key + ".zip", origin=details["url"],
                                  file_hash=details["hash"], cache_subdir=target, extract=True))
        candidates = [extracted, extracted / key, extracted.parent / key]
        folder = next((path for path in candidates if (path / "config.json").is_file()), None)
        if folder is None:
            raise RuntimeError(f"StarDist model config not found in {extracted}") from error
        model = model_class(config=None, name=folder.name, basedir=folder.parent)
        model.basedir = None
        return model


def segment_cells(
    image_path: str,
    output_path: str = "./nuclei_labels.npz",
    polys_path: str = None,
    plot_path: str = None,
    model_type: str = "2D_versatile_he",
    tile_size: int = 4096,
    overlap: int = 128,
    context: int = 128,
    min_area: int = 20,
    prob_thresh: float = None,
    nms_thresh: float = None,
    normalize_pct: tuple = (5, 95),
) -> str:
    """Segment nuclei from an H&E image using StarDist (10x official style).

    Mirrors 10x NucleiSegmentation.py: percentile normalization (5,95),
    StarDist 2D_versatile_he, predict_instances_big with
    prob_thresh/nms_thresh low, saves BOTH the label array (NPZ) and the
    polygon dict (numeric NPZ, key 'coord') used by :func:`bin_to_cell` for
    geopandas sjoin 'within'.

    Args:
        image_path: Path to the H&E image (btf/tif/png/jpg).
        output_path: Where to save the label array (NPZ with key 'labels').
        polys_path: Where to save the polygon dict (numeric NPZ). Defaults to
            output_path with '_polys.npz' suffix.
        model_type: StarDist model name ('2D_versatile_he' or '2D_versatile_fluo').
        tile_size: Tile size for big-image inference.
        overlap: Tile overlap in pixels.
        context: Context in pixels for big-image inference (official 128).
        min_area: Minimum nucleus area (px) to keep.
        prob_thresh: Override StarDist probability threshold (None = official 0.01).
        nms_thresh: Override NMS threshold (None = official 0.001).
        normalize_pct: Percentile normalization (official 5,95).

    Returns:
        str: Log with nucleus count and output paths.
    """
    import os as _os
    # GPU 优先: 默认用 GPU, 但 cuDNN 与驱动不兼容的机器会在编译卷积时
    # 抛 Autotuner error (进程级污染, 事后无法回退), 所以必须在 TF 首次
    # 初始化之前探测。STARDIST_DEVICE=cpu 可强制走 CPU。
    _force_cpu = _os.environ.get("STARDIST_DEVICE", "").lower() == "cpu"
    if _force_cpu or not _probe_gpu():
        _os.environ["CUDA_VISIBLE_DEVICES"] = ""
    # 模型已下载到本地缓存后, 跳过联网检查 (避免网络抖动导致失败)
    if _os.path.exists(_os.path.expanduser("~/.keras/models/StarDist2D")):
        _os.environ.setdefault("STARDIST_OFFLINE", "1")
        _os.environ.setdefault("CSBDEEP_OFFLINE", "1")

    try:
        from stardist.models import StarDist2D
    except ImportError:
        return ("Error: stardist not installed. Install with: "
                "pip install stardist tensorflow")

    import numpy as np

    if polys_path is None:
        polys_path = _os.path.splitext(output_path)[0] + "_polys.npz"

    # 按扩展名选择读取器: tif/tiff 用 tifffile, 其他(png/jpg)用 imageio
    ext = _os.path.splitext(image_path)[1].lower()
    if ext in (".tif", ".tiff", ".btf"):
        from tifffile import imread
        img = imread(image_path)
    else:
        import imageio.v2 as imageio
        img = imageio.imread(image_path)
    # 灰度图转 RGB 或保持; StarDist 需要 (H, W, C)
    if img.ndim == 2:
        img = np.stack([img] * 3, axis=-1)
    elif img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]
    # 官方 normalize (percentile 5-95), 与 10x NucleiSegmentation.py 一致
    from csbdeep.utils import normalize
    img = normalize(img, normalize_pct[0], normalize_pct[1])

    model = _load_stardist_pretrained(StarDist2D, model_type)
    # 官方阈值: prob_thresh=0.01, nms_thresh=0.001
    if prob_thresh is None:
        prob_thresh = 0.01
    if nms_thresh is None:
        nms_thresh = 0.001

    # stardist >= 0.9: predict_instances(_big) 是模型方法, 需显式 axes
    # 小图(<=tile_size)直接用 predict_instances, 避免 big 路径的 context 开销丢核
    h_img, w_img = img.shape[:2]
    if h_img <= tile_size and w_img <= tile_size:
        labels, polys = model.predict_instances(
            img, axes="YXC", prob_thresh=prob_thresh, nms_thresh=nms_thresh
        )
    else:
        labels, polys = model.predict_instances_big(
            img, axes="YXC", block_size=tile_size, min_overlap=overlap,
            context=context, prob_thresh=prob_thresh, nms_thresh=nms_thresh,
            normalizer=None, n_tiles=(4, 4, 1),
        )

    labels, polys = filter_nuclei_instances(labels, polys, min_area)
    n_nuclei = len(polys["coord"])
    _os.makedirs(_os.path.dirname(output_path) or ".", exist_ok=True)
    save_nuclei_polygons(polys_path, polys)
    np.savez_compressed(output_path, labels=labels)

    # 核分割检查图 (私有绘图函数, 函数内保存 PNG + 写进日志)
    if plot_path is None:
        plot_path = os.path.splitext(output_path)[0] + "_crop.png"
    try:
        _plot_nuclei_crop(labels, image_path, plot_path,
                          title=f"{n_nuclei} nuclei")
        crop_note = f"; crop plot -> {plot_path}"
    except Exception as exc:
        crop_note = f"; crop plot FAILED: {exc}"
    return (f"Segmented {n_nuclei} nuclei from {_os.path.basename(image_path)} "
            f"-> labels {output_path}, polys {polys_path}{crop_note}")


def bin_to_cell(
    adata_path: str,
    polys_path: str,
    output_path: str = "./cells.h5ad",
    barcode_map_path: str = None,
    image_path: str = None,
    plot_path: str = None,
    expand_px: float = 29.0,
    allow_legacy_pickle: bool = False,
) -> str:
    """Assign 2 um bins to nuclei via geopandas sjoin (10x official + expansion).

    Mirrors 10x NucleiSegmentation.py: build Polygon GeoDataFrame from StarDist
    polys, build Point GeoDataFrame from bin coordinates, sjoin bins to nuclei.

    **Nucleus expansion (10x Space Ranger default 8 um):** the raw nucleus
    polygon covers only the nucleus; without expansion each cell gets only the
    transcripts within the nucleus (median ~3 bins, very low UMI). 10x Space
    Ranger expands nuclei by ``--nucleus-expansion-distance-micron`` (default
    8 um); bin2cell/ENACT also expand by N bins. Here we buffer the polygons by
    ``expand_px`` (pixels) before sjoin, and assign overlapping bins to the
    nearest nucleus.

    Args:
        adata_path: Path to the loaded Visium HD h5ad (bins x genes).
        polys_path: Path to numeric NPZ with StarDist polys (key 'coord').
        output_path: Where to save the cell-level AnnData.
        barcode_map_path: Optional path for the official Nuclei_Barcode_Map.csv.
        image_path: Optional H&E image for overlay plot.
        plot_path: Optional overlay plot path.
        expand_px: Polygon buffer distance in fullres pixels (10x default 8 um
            ~= 29 px at 0.2738 um/px). 0 disables expansion (strict within).

    Returns:
        str: Log with cell count and output path.
    """
    import numpy as np
    import anndata as ad
    import pandas as pd
    from scipy import sparse
    import geopandas as gpd
    from shapely.geometry import Polygon, Point

    adata = ad.read_h5ad(adata_path)
    polys = load_nuclei_polygons(polys_path, allow_legacy_pickle=allow_legacy_pickle)
    if not np.isfinite(expand_px) or expand_px < 0:
        raise ValueError("expand_px must be finite and nonnegative")
    if not adata.obs_names.is_unique:
        raise ValueError("Bin barcodes must be unique")

    if "spatial" not in adata.obsm:
        return "Error: adata lacks obsm['spatial']; cannot assign bins to cells"

    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    if coords.shape != (adata.n_obs, 2) or not np.isfinite(coords).all():
        raise ValueError("Bin spatial coordinates must be finite (n_bins, 2) values")
    geometries = [Polygon(vertices[::-1].T) for vertices in polys["coord"]]
    if not geometries or any(g.is_empty or not g.is_valid or g.area <= 0 for g in geometries):
        raise ValueError("Provide at least one valid nonempty nucleus polygon")
    gdf = gpd.GeoDataFrame(geometry=geometries)
    gdf["object_id"] = np.arange(1, len(gdf) + 1)
    gdf["id"] = [f"ID_{i}" for i in gdf["object_id"]]
    gdf["nucleus_x"] = gdf.geometry.centroid.x
    gdf["nucleus_y"] = gdf.geometry.centroid.y
    if expand_px:
        gdf["geometry"] = gdf.geometry.buffer(expand_px)
    bins = gpd.GeoDataFrame({"bin_row": np.arange(adata.n_obs), "index": adata.obs_names},
                           geometry=gpd.points_from_xy(coords[:, 0], coords[:, 1]))
    candidates = gpd.sjoin(bins, gdf, how="inner", predicate="intersects")
    if candidates.empty:
        raise ValueError("No bins overlap the nucleus polygons; check image and coordinate alignment")
    candidates["is_overlap"] = candidates.duplicated("bin_row", keep=False)
    candidates["_distance2"] = ((candidates.geometry.x - candidates["nucleus_x"]) ** 2 +
                                (candidates.geometry.y - candidates["nucleus_y"]) ** 2)
    result = (candidates.sort_values(["bin_row", "_distance2", "object_id"], kind="stable")
              .drop_duplicates("bin_row", keep="first").copy())
    result["is_within_polygon"] = True
    if barcode_map_path:
        os.makedirs(os.path.dirname(barcode_map_path) or ".", exist_ok=True)
        result[["index", "id", "is_within_polygon", "is_overlap"]].to_csv(barcode_map_path, index=False)

    rows_idx = result["bin_row"].to_numpy(dtype=np.int64)
    cell_ids = np.sort(result["object_id"].unique()).astype(np.int64)
    cols_idx = np.searchsorted(cell_ids, result["object_id"].to_numpy())
    assignment = sparse.csr_matrix((np.ones(len(rows_idx), dtype=np.int64), (rows_idx, cols_idx)),
                                   shape=(adata.n_obs, len(cell_ids)))
    cell_gene = assignment.T @ adata.X
    counts = np.bincount(cols_idx, minlength=len(cell_ids))
    # Each assigned bin occurs once, so every gene's assigned count is conserved.
    expected = np.asarray(adata.X[rows_idx].sum(axis=0)).ravel()
    actual = np.asarray(cell_gene.sum(axis=0)).ravel()
    if not np.allclose(actual, expected):
        raise RuntimeError("Bin aggregation failed per-gene count conservation")
    cell_adata = ad.AnnData(X=cell_gene,
                           obs=pd.DataFrame(index=[f"cell_{cid}" for cid in cell_ids]),
                           var=adata.var.copy())
    cell_adata.obsm["spatial"] = np.asarray(assignment.T @ coords) / counts[:, None]
    cell_adata.obs["object_id"] = cell_ids
    cell_adata.obs["bin_count"] = counts
    cell_adata.uns["bin_to_cell"] = {
        "expand_px": float(expand_px), "n_bins_total": int(adata.n_obs),
        "n_bins_assigned": int(len(result)), "n_bins_unassigned": int(adata.n_obs - len(result)),
        "n_cells": int(len(cell_ids)), "frac_bins_assigned": len(result) / adata.n_obs,
        "overlap_rule": "nearest_original_nucleus_centroid_then_object_id",
    }
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    cell_adata.write(output_path, compression="lzf")

    # 细胞 overlay 图 (HE + 质心, 私有绘图函数, 函数内保存 PNG + 写进日志)
    if plot_path is None:
        plot_path = os.path.splitext(output_path)[0] + "_overlay.png"
    if image_path:
        try:
            _plot_cells_overlay(cell_adata, image_path, plot_path,
                                title=f"{len(cell_ids)} cells (StarDist, expand {expand_px}px)",
                                geometry_path=polys_path, geometry_kind="stardist",
                                edge_color="tab:blue", linewidth=0.12, allow_legacy_pickle=allow_legacy_pickle)
            overlay_note = f"; overlay -> {plot_path}"
        except Exception as exc:
            overlay_note = f"; overlay FAILED: {exc}"
    else:
        overlay_note = ""
    return (f"Aggregated {adata.n_obs} bins -> {len(cell_ids)} cells "
            f"({len(result)} bins assigned, expand_px={expand_px}) -> {output_path}{overlay_note}")


def run_proseg(
    visiumhd_dir: str,
    labels_npz: str,
    output_dir: str = "./proseg_out",
    voxel_size: float = None,
    burnin_voxel_size: float = None,
    samples: int = None,
    burnin_samples: int = None,
    recorded_samples: int = None,
    no_diffusion: bool = None,
    cellpose_scale: float = 0.2738,
    quad_size: int = 150,
    max_transcript_nucleus_distance: int = 60,
    nthreads: int = 32,
    proseg_executable: str = "proseg",
) -> str:
    """Run Proseg Bayesian cell segmentation on Visium HD binned data.

    Self-contained: takes the NPZ from :func:`segment_cells`, converts the
    labels to the uint32 NPY format Proseg requires internally, reads the 2 um
    binned transcripts, and runs MCMC over voxels. Outputs a cell-level counts
    matrix (cell x gene) + metadata + spatialdata zarr.

    Args:
        visiumhd_dir: Path to Space Ranger output dir (contains binned_outputs/).
        labels_npz: Path to NPZ with 'labels' from :func:`segment_cells`.
        output_dir: Where to write proseg outputs (zarr, counts, meta).
        voxel_size: Final x/y voxel size in microns (None = auto by tissue size).
        burnin_voxel_size: Burn-in voxel size (None = 2x voxel_size).
        samples: Number of post-burnin samples (None = auto by transcript count).
        burnin_samples: Burn-in samples (None = same as samples).
        recorded_samples: Samples recorded for output (None = samples//2).
        no_diffusion: Disable transcript diffusion between cells (None = auto by memory).
        cellpose_scale: Microns per pixel of the mask image (Visium HD fullres ~0.2738).
        quad_size: Proseg quad-tree split threshold (default 150).
        max_transcript_nucleus_distance: Max distance transcripts assigned to nuclei.
        nthreads: CPU threads for Proseg.

    Returns:
        str: Log with cell count and output paths.
    """
    import os as _os
    import subprocess
    import time
    import numpy as np

    _os.makedirs(output_dir, exist_ok=True)

    # --- 自包含: NPZ(labels) -> uint32 NPY (Proseg 强制要求) ---
    masks_path = _os.path.join(output_dir, "proseg_masks.npy")
    if not np.isfinite(cellpose_scale) or cellpose_scale <= 0:
        raise ValueError("cellpose_scale must be finite and positive")
    with np.load(labels_npz, allow_pickle=False) as archive:
        labels = archive["labels"]
    if labels.ndim != 2 or not np.issubdtype(labels.dtype, np.integer) or (labels < 0).any():
        raise ValueError("Proseg masks must be nonnegative integer labels in two dimensions")
    if labels.max(initial=0) > np.iinfo(np.uint32).max:
        raise ValueError("Proseg mask IDs exceed uint32")
    if labels.dtype != np.uint32:
        labels = labels.astype(np.uint32)
    np.save(masks_path, labels)

    # --- 自动推断参数 (不写死, 换样本自适应) ---
    # Proseg read_visium_data 期望 path 指向包含 square_002um 的目录 (= binned_outputs)
    binned_dir = _os.path.join(visiumhd_dir, "binned_outputs")
    if not _os.path.isdir(binned_dir):
        binned_dir = visiumhd_dir  # 容错: 如果直接传了 binned_outputs
    import h5py
    h5 = _os.path.join(binned_dir, "square_002um",
                       "filtered_feature_bc_matrix.h5")
    transcripts = 0
    img_w, img_h = 0, 0
    if _os.path.exists(h5):
        with h5py.File(h5, "r") as f:
            grp = f["matrix"]
            dataset = grp["data"]
            transcripts = sum(int(dataset[i:i + 1_000_000].sum(dtype=np.int64))
                              for i in range(0, len(dataset), 1_000_000))
    else:
        raise FileNotFoundError(h5)
    img_path = None
    for cand in (_os.path.join(visiumhd_dir, "tissue_image.btf"),
                 _os.path.join(visiumhd_dir, "tissue_image.tiff"),
                 _os.path.join(visiumhd_dir, "spatial", "cytassist_image.tiff")):
        if _os.path.exists(cand):
            img_path = cand
            break
    if img_path:
        import tifffile
        with tifffile.TiffFile(img_path) as tf:
            page = tf.pages[0]
            img_w, img_h = page.shape[1], page.shape[0]
    mem_gb = 32.0
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable"):
                    mem_gb = int(line.split()[1]) / 1e6
                    break
    except Exception:
        pass

    tissue_size_px = max(img_w, img_h)
    if voxel_size is None:
        voxel_size = 1.0 if tissue_size_px < 6000 else 2.0
    if burnin_voxel_size is None:
        burnin_voxel_size = voxel_size * 2
    if samples is None:
        samples = 200 if transcripts < 5e7 else (150 if transcripts < 1e8 else 100)
    if burnin_samples is None:
        burnin_samples = samples
    if recorded_samples is None:
        recorded_samples = max(1, min(samples, samples // 2))

    # 内存估算: voxel 数 * 状态 + 转录本状态
    mpp = cellpose_scale
    w_um, h_um = img_w * mpp, img_h * mpp
    n_voxels = (w_um / voxel_size) * (h_um / voxel_size) * 4  # z 层 4
    est_gb = (n_voxels * 60 + transcripts * 100 + transcripts * 4) / 1e9
    if no_diffusion is None:
        no_diffusion = est_gb > mem_gb * 0.6

    # --- 组装 proseg 命令 ---
    cmd = [str(proseg_executable), "--visiumhd",
           "--cellpose-masks", masks_path,
           "--cellpose-scale", str(cellpose_scale),
           "--voxel-size", str(voxel_size),
           "--burnin-voxel-size", str(burnin_voxel_size),
           "--samples", str(samples),
           "--burnin-samples", str(burnin_samples),
           "--recorded-samples", str(recorded_samples),
           "--quad-size", str(quad_size),
           "--max-transcript-nucleus-distance", str(max_transcript_nucleus_distance),
           "--nthreads", str(nthreads),
           "--output-spatialdata", _os.path.join(output_dir, "proseg_out.zarr"),
           "--output-counts", _os.path.join(output_dir, "proseg_counts.mtx.gz"),
           "--output-cell-metadata", _os.path.join(output_dir, "proseg_meta.csv.gz"),
           "--output-gene-metadata", _os.path.join(output_dir, "proseg_genes.csv.gz"),
           binned_dir]
    if no_diffusion:
        cmd.append("--no-diffusion")

    env = _os.environ.copy()
    env["PATH"] = _os.path.expanduser("~/.cargo/bin") + _os.pathsep + env.get("PATH", "")
    t0 = time.time()
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=7200)
    elapsed = round(time.time() - t0, 1)
    if proc.returncode != 0:
        return (f"Error: proseg failed rc={proc.returncode} ({elapsed}s)\n"
                f"stderr tail: {proc.stderr[-2000:]}")

    # --- 读取结果 -> cell h5ad ---
    import numpy as np
    import pandas as pd
    import anndata as ad
    from scipy.io import mmread
    from scipy.sparse import csr_matrix

    counts_path = _os.path.join(output_dir, "proseg_counts.mtx.gz")
    meta_path = _os.path.join(output_dir, "proseg_meta.csv.gz")
    genes_path = _os.path.join(output_dir, "proseg_genes.csv.gz")
    cell_adata = read_proseg_counts(counts_path, meta_path, genes_path,
                                   microns_per_pixel=cellpose_scale)
    M = cell_adata.X

    out_h5ad = _os.path.join(output_dir, "proseg_cells.h5ad")
    cell_adata.write(out_h5ad, compression="lzf")

    # Proseg 细胞 overlay 图 (HE + 质心)
    overlay_note = ""
    if img_path:
        try:
            plot_p = _os.path.join(output_dir, "proseg_cells_overlay.png")
            geom_p = _os.path.join(output_dir, "proseg_out.zarr",
                                   "shapes", "cell_boundaries", "shapes.parquet")
            _plot_cells_overlay(cell_adata, img_path, plot_p,
                                title=f"{M.shape[0]} cells (Proseg)",
                                geometry_path=geom_p, geometry_kind="proseg",
                                edge_color="tab:red", linewidth=0.12)
            overlay_note = f"; overlay -> {plot_p}"
        except Exception as exc:
            overlay_note = f"; overlay FAILED: {exc}"

    return (f"Proseg done: {M.shape[0]} cells x {M.shape[1]} genes ({elapsed}s, "
            f"voxel={voxel_size}um samples={samples} no_diffusion={no_diffusion}, "
            f"est_mem={est_gb:.1f}G avail={mem_gb:.0f}G) -> {out_h5ad}{overlay_note}")


def load_xenium(
    matrix_path: str,
    coordinates_path: str = None,
    sample_id: str = None,
    output_path: str = "./xenium_cells.h5ad",
) -> str:
    """Load Xenium cell-level matrix (platform segmentation is ground truth).

    Xenium outputs a cell x gene matrix with platform-provided segmentation.
    No re-segmentation or deconvolution is needed.

    Args:
        matrix_path: Path to cell_feature_matrix.h5 or cell-level h5ad.
        coordinates_path: Optional CSV with cell_id,x,y (if h5ad lacks obsm['spatial']).
        sample_id: Optional sample identifier.
        output_path: Where to save the cell-level AnnData.

    Returns:
        str: Log.
    """
    import scanpy as sc
    import anndata as ad

    steps = []
    if matrix_path.endswith(".h5ad"):
        adata = ad.read_h5ad(matrix_path)
        steps.append(f"Loaded h5ad: {adata.shape[0]} cells x {adata.shape[1]} genes")
    else:
        adata = sc.read_10x_h5(matrix_path)
        adata.var_names_make_unique()
        steps.append(f"Loaded 10x h5: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Xenium 坐标通常在 obs 的 pxl_col/pxl_row 或 obsm['spatial']
    if "spatial" not in adata.obsm:
        for xk, yk in (("pxl_col_in_fullres", "pxl_row_in_fullres"),
                       ("x", "y"), ("x_centroid", "y_centroid")):
            if xk in adata.obs and yk in adata.obs:
                adata.obsm["spatial"] = adata.obs[[xk, yk]].values.astype(float)
                steps.append(f"Recovered spatial from obs {xk}/{yk}")
                break

    if coordinates_path and os.path.exists(coordinates_path) and "spatial" not in adata.obsm:
        import pandas as pd
        pos = pd.read_csv(coordinates_path)
        col_map = {}
        for i, c in enumerate(pos.columns):
            cl = str(c).strip().lower()
            if cl in ("barcode", "bc", "cell_id", "id"):
                col_map["barcode"] = c
            elif cl in ("pxl_col_in_fullres", "pxl_col", "x_pixel", "x", "x_centroid"):
                col_map["x"] = c
            elif cl in ("pxl_row_in_fullres", "pxl_row", "y_pixel", "y", "y_centroid"):
                col_map["y"] = c
        if "x" in col_map and "y" in col_map:
            if "barcode" in col_map:
                pos.index = pos[col_map["barcode"]].astype(str).values
                adata.obsm["spatial"] = pos.loc[adata.obs_names, [col_map["x"], col_map["y"]]].astype(float).values
            elif len(pos) == adata.n_obs:
                adata.obsm["spatial"] = pos[[col_map["x"], col_map["y"]]].astype(float).values
            steps.append(f"Attached coordinates from {os.path.basename(coordinates_path)}")

    if sample_id:
        adata.obs["sample_id"] = sample_id

    adata.uns["spatial_meta"] = {
        "source": matrix_path,
        "platform": "xenium",
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "has_coords": "spatial" in adata.obsm,
    }

    adata.write(output_path, compression="lzf")
    steps.append(f"Saved Xenium cells to {output_path}")
    return "\n".join(steps)


def load_atera(
    matrix_path: str,
    sample_id: str = None,
    output_path: str = "./atera_cells.h5ad",
    *,
    roi: tuple = None,
    use_official_filter: bool = True,
    cell_type_key: str = "azimuth_coarse",
    chunk_rows: int = 2048,
) -> str:
    """Load official Atera AnnData Zarr or a previously prepared h5ad.

    Atera is 10x's whole-transcriptome in situ platform. Official output is
    cell-level AnnData + Zarr segmentation. No re-segmentation or deconvolution.

    Args:
        matrix_path: Output directory, experiment.spatial, cell_feature_matrix.zarr.zip, or h5ad.
        sample_id: Optional sample identifier.
        output_path: Where to save the cell-level AnnData.
        roi: For native Zarr, (xmin, ymin, xmax, ymax), in micrometres.
        use_official_filter: For native Zarr, retain obs['filtered'] == True.
        cell_type_key: Official annotation field to copy to cell_type; None keeps labels unset.
        chunk_rows: Rows per read when converting the native sparse matrix.

    Returns:
        str: Log.
    """
    import anndata as ad

    matrix_path = str(matrix_path)
    if not matrix_path.endswith(".h5ad"):
        try:
            from .atera import load_atera_zarr
        except ImportError:
            from atera import load_atera_zarr
        return load_atera_zarr(matrix_path, output_path, sample_id=sample_id, roi=roi,
                               use_official_filter=use_official_filter,
                               cell_type_key=cell_type_key, chunk_rows=chunk_rows)
    if roi is not None:
        raise ValueError("Pass the native Atera Zarr for ROI extraction in micrometres")
    steps = []
    adata = ad.read_h5ad(matrix_path)
    steps.append(f"Loaded Atera AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Atera 坐标通常在 obsm['spatial']; 若无则从 obs 恢复
    if "spatial" not in adata.obsm:
        for xk, yk in (("x", "y"), ("x_centroid", "y_centroid"),
                       ("pxl_col_in_fullres", "pxl_row_in_fullres")):
            if xk in adata.obs and yk in adata.obs:
                adata.obsm["spatial"] = adata.obs[[xk, yk]].values.astype(float)
                steps.append(f"Recovered spatial from obs {xk}/{yk}")
                break

    if sample_id:
        adata.obs["sample_id"] = sample_id

    meta = dict(adata.uns.get("spatial_meta", {}))
    if cell_type_key is not None and cell_type_key in adata.obs:
        adata.obs["cell_type"] = adata.obs[cell_type_key].astype("category")
        meta["cell_type_source"] = cell_type_key
    meta.update({
        "source": matrix_path,
        "platform": "atera",
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "has_coords": "spatial" in adata.obsm,
    })
    adata.uns["spatial_meta"] = meta

    adata.write(output_path, compression="lzf")
    steps.append(f"Saved Atera cells to {output_path}")
    return "\n".join(steps)
