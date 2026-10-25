"""

Spatial transcriptomics (10x Visium) analysis tools.



Implements a standard Visium analysis workflow with per-step visualization:

  Stage 1: data loading, QC filtering, normalization

  Stage 2: clustering and spatial domain detection

  Stage 3: spatially variable gene detection

  Stage 4: cell-type deconvolution (SPOTlight / DestVI)

  Stage 5: cell-cell communication inference (CellChat v2)

  Stage 6: neighborhood enrichment analysis



Every tool saves its own plots (PNG) under the output directory and returns

the plot paths in the step log.

"""



import os



import numpy as np

import pandas as pd

import scanpy as sc






def _run_r_script(executable, script, output_dir):
    """Keep R diagnostics on disk while the subprocess is running."""
    import subprocess
    from pathlib import Path
    stdout_path = Path(output_dir) / "r_stdout.log"
    stderr_path = Path(output_dir) / "r_stderr.log"
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        result = subprocess.run([executable, script], stdout=stdout, stderr=stderr, timeout=3600)
    result.stdout = stdout_path.read_text(encoding="utf-8", errors="replace")
    result.stderr = stderr_path.read_text(encoding="utf-8", errors="replace")
    return result


def _r_literal(value):
    """Encode a path or label as a quoted R string, including Windows paths."""
    import json
    return json.dumps(str(value), ensure_ascii=False)


def _raw_count_adata(adata, source):
    """Choose preserved count data and validate the count-model input."""
    from scipy.sparse import issparse
    candidates = []
    if "counts" in adata.layers:
        candidate = adata.copy()
        candidate.X = adata.layers["counts"].copy()
        candidates.append(candidate)
    if adata.raw is not None:
        candidate = adata.raw.to_adata()
        candidate.uns.pop("log1p", None)
        candidates.append(candidate)
    if "log1p" not in adata.uns:
        candidates.append(adata.copy())
    for candidate in candidates:
        values = candidate.X.data if issparse(candidate.X) else np.asarray(candidate.X).ravel()
        if np.isfinite(values).all() and (values >= 0).all() and np.allclose(values, np.rint(values), rtol=0, atol=1e-6):
            return candidate
    raise ValueError(f"{source}: provide nonnegative integer raw counts in layers['counts'], raw.X, or X")


def _plot_setup(output_dir: str) -> None:

    """Ensure matplotlib uses the Agg backend and the output dir exists."""

    import matplotlib



    matplotlib.use("Agg")

    os.makedirs(output_dir, exist_ok=True)

    # FIX 2026-08-26: 固定字体规范 (用户要求) — 英文 Arial, 中文 SimHei.
    # 优先从常见路径 addfont, 缺失时静默跳过 (fallback 默认字体).
    import matplotlib.font_manager as _fm
    import matplotlib.pyplot as _plt
    for _p in ("/usr/share/fonts/truetype/custom/arial.ttf",
               "/usr/share/fonts/truetype/custom/arialbd.ttf",
               "/usr/share/fonts/truetype/custom/simhei.ttf",
               "C:\\Windows\\Fonts\\arial.ttf",
               "C:\\Windows\\Fonts\\simhei.ttf"):
        try:
            if os.path.exists(_p):
                _fm.fontManager.addfont(_p)
        except Exception:
            pass
    _plt.rcParams["font.family"] = ["Arial", "SimHei"]
    _plt.rcParams["axes.unicode_minus"] = False





def _categorical_cmap(n_cats):
    """Return a list of n distinct colors (no repeats up to 60, then HSV cycle).

    tab20/tab20b/tab20c each have 20 colors -> 60 distinct before cycling.
    For >60 categories, fall back to HSV with golden-ratio spacing (visually distinct).
    """
    import matplotlib.pyplot as plt
    if n_cats <= 60:
        cols = []
        for name in ("tab20", "tab20b", "tab20c"):
            cols.extend(plt.get_cmap(name).colors)
        return cols[:n_cats]
    # >60: HSV golden-ratio spacing
    import colorsys
    phi = 0.618033988749895
    cols = []
    for i in range(n_cats):
        h = (i * phi) % 1.0
        r, g, b = colorsys.hsv_to_rgb(h, 0.85, 0.85)
        cols.append((r, g, b))
    return cols


def detect_species(adata, n_marker=6):
    """Detect species from gene symbols in an AnnData object.

    Returns 'mouse' or 'human'. Signals (priority):
      1. Mitochondrial gene prefix: mt- (mouse, lowercase) vs MT- (human, uppercase)
      2. Known marker gene hit rate (mouse vs human symbols)
      3. TitleCase fraction: >80% of genes TitleCase (Xkr4-like) -> mouse;
         all-uppercase symbols (XKR4-like) -> human
    If ambiguous, returns 'unknown' (LLM should decide).
    """
    import numpy as np
    genes = list(adata.var_names)
    if not genes:
        return "unknown"

    # Signal 1: mitochondrial prefix
    mt_mouse = sum(1 for g in genes if g.startswith("mt-"))
    mt_human = sum(1 for g in genes if g.startswith("MT-"))
    if mt_mouse > 0 and mt_human == 0:
        return "mouse"
    if mt_human > 0 and mt_mouse == 0:
        return "human"

    # Signal 2: known markers (mouse symbol style vs human)
    markers = {
        "mouse": ["Tgfbr1", "Bmp2", "Cdh1", "App", "Ptprc", "Cd3d", "Lyz1", "Defa21"],
        "human": ["TGFBR1", "BMP2", "CDH1", "APP", "PTPRC", "CD3D", "LYZ", "DEFA5"],
    }
    gene_set = set(genes)
    hits = {sp: sum(1 for m in ms if m in gene_set) for sp, ms in markers.items()}
    if hits["mouse"] >= 2 and hits["human"] == 0:
        return "mouse"
    if hits["human"] >= 2 and hits["mouse"] == 0:
        return "human"

    # Signal 3: TitleCase fraction
    title = sum(1 for g in genes if g[:1].isupper() and not g.isupper())
    upper = sum(1 for g in genes if g.isupper())
    if title > 0 and upper > 0:
        frac = title / (title + upper)
        if frac > 0.8:
            return "mouse"
        if frac < 0.2:
            return "human"
    return "unknown"


def _spatial_scatter(adata, color, ax, title, cmap="viridis"):

    """Plot spots colored by a value/annotation on the spatial coordinates."""

    import matplotlib.pyplot as plt



    if "spatial" not in adata.obsm:

        return

    coords = adata.obsm["spatial"]

    x, y = coords[:, 0], coords[:, 1]

    if color in adata.obs:

        vals = adata.obs[color]

        if pd.api.types.is_categorical_dtype(vals) or vals.dtype == object:

            cats = vals.cat.categories if pd.api.types.is_categorical_dtype(vals) else pd.unique(vals)

            # FIX 2026-08-26: 分类配色不重复 (tab20 只有 20 色, >20 类会循环重复)
            cmap_colors = _categorical_cmap(len(cats))

            color_map = {c: cmap_colors[i] for i, c in enumerate(cats)}

            colors = [color_map.get(v, "grey") for v in vals]

            scatter = ax.scatter(x, y, c=colors, s=8, alpha=0.8)

            handles = [plt.Line2D([0], [0], marker="o", color="w",
                                  markerfacecolor=color_map[c], markersize=6,
                                  label=str(c)) for c in cats]

            ax.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc="upper left",
                      fontsize=7, frameon=False, title=color, title_fontsize=8)

        else:

            scatter = ax.scatter(x, y, c=vals, s=8, alpha=0.8, cmap=cmap)

            plt.colorbar(scatter, ax=ax, fraction=0.046)

    elif color in adata.var_names:

        values = adata[:, color].X
        values = values.toarray() if hasattr(values, "toarray") else np.asarray(values)
        scatter = ax.scatter(x, y, c=values.ravel(), s=8, alpha=0.8, cmap=cmap)

        plt.colorbar(scatter, ax=ax, fraction=0.046)

    else:

        ax.scatter(x, y, s=8, alpha=0.8)

    ax.set_title(title)

    ax.set_xlabel("x")

    ax.set_ylabel("y")

    ax.set_aspect("equal")





# ============================================================================

# Stage 1: Data loading, QC and normalization

# ============================================================================





def load_visium_data(

    counts_file: str,

    coordinates_file: str = None,

    sample_id: str = None,

    output_path: str = "./visium_loaded.h5ad",

    plot_path: str = None,

) -> str:

    """Load Visium data into a standardized AnnData object.



    This tool does NOT guess file formats. The caller (LLM agent) is

    responsible for identifying the input data format first and passing

    explicit file paths. Supported inputs:



      1. h5ad file: counts_file = path to .h5ad (must contain raw counts in X,

         optionally obsm['spatial']).

      2. 10x matrix directory: counts_file = directory containing

         matrix.mtx(.gz) + features.tsv(.gz) + barcodes.tsv(.gz).

      3. Wide counts CSV: counts_file = path to CSV with rows=barcodes,

         cols=genes (e.g. gastric GSE counts).



    If the data is not one of these formats, do NOT call this tool. Instead

    inspect the files yourself (use run_python_repl), convert them to a

    supported format (or to h5ad), then call this tool.



    Args:

        counts_file: Path to h5ad, 10x matrix directory, or wide counts CSV.

        coordinates_file: Optional path to coordinate CSV (barcode + x/y or

            Space Ranger tissue_positions format). If the h5ad already has

            obsm['spatial'] this is not needed.

        sample_id: Optional sample identifier stored in obs['sample_id'].

        output_path: Path where the loaded AnnData will be saved.

        plot_path: Path where the spot-layout plot will be saved.



    Returns:

        str: Step-by-step log describing the load result.

    """

    import scanpy as sc



    steps = []

    if not os.path.exists(counts_file):

        return f"Error: counts_file not found: {counts_file}"



    # ---- Load counts ----

    if counts_file.endswith(".h5ad"):

        adata = sc.read_h5ad(counts_file)

        steps.append(f"Loaded h5ad: {adata.shape[0]} spots x {adata.shape[1]} genes")

        # Try to recover coordinates from obs if not in obsm

        if "spatial" not in adata.obsm:

            for xk, yk in (("pixel_x", "pixel_y"), ("pxl_col_in_fullres", "pxl_row_in_fullres"),

                           ("x", "y"), ("array_x", "array_y")):

                if xk in adata.obs and yk in adata.obs:

                    adata.obsm["spatial"] = adata.obs[[xk, yk]].values.astype(float)

                    steps.append(f"Recovered coordinates from obs columns {xk}/{yk}")

                    break

    elif os.path.isdir(counts_file):

        # 10x matrix directory

        import scipy.io as sio



        mtx = os.path.join(counts_file, "matrix.mtx")

        feat = os.path.join(counts_file, "features.tsv")

        bar = os.path.join(counts_file, "barcodes.tsv")

        paths = [mtx, feat, bar]
        for i, path in enumerate(paths):
            if not os.path.exists(path):
                if os.path.exists(path + ".gz"):
                    paths[i] = path + ".gz"
                else:
                    raise FileNotFoundError(path)
        mtx, feat, bar = paths
        import gzip
        with (gzip.open(mtx, "rb") if mtx.endswith(".gz") else open(mtx, "rb")) as handle:
            counts = sio.mmread(handle).T.tocsr()
        ft = pd.read_csv(feat, sep="\t", header=None)
        bc = pd.read_csv(bar, sep="\t", header=None)
        adata = sc.AnnData(X=counts)
        adata.var_names = ft.iloc[:, 1 if ft.shape[1] >= 2 else 0].astype(str).values
        adata.var_names_make_unique()
        adata.obs_names = bc.iloc[:, 0].astype(str).values
        adata.obs_names_make_unique()

        steps.append(f"Loaded 10x matrix: {adata.shape[0]} spots x {adata.shape[1]} genes")

    else:

        # Wide counts CSV (rows=barcodes, cols=genes)

        import gzip



        opener = gzip.open if counts_file.endswith(".gz") else open

        with opener(counts_file, "rt", errors="replace") as f:

            df = pd.read_csv(f, index_col=0)

        df.index = df.index.astype(str)

        adata = sc.AnnData(X=df.values.astype(float))

        adata.obs_names = df.index.values

        adata.var_names = df.columns.astype(str).values

        adata.var_names_make_unique()

        steps.append(f"Loaded wide counts CSV: {adata.shape[0]} spots x {adata.shape[1]} genes")



    # ---- Standardize ----

    if "gene_symbols" not in adata.var:

        adata.var["gene_symbols"] = adata.var_names.astype(str).values

    adata.obs["barcode"] = adata.obs_names.astype(str).values



    # ---- Attach coordinates if provided ----

    if coordinates_file and os.path.exists(coordinates_file):
        raw = pd.read_csv(coordinates_file, header=None)
        # --- Detect header by column-name semantics (not by column count) ---
        first_row = [str(v).strip().lower() for v in raw.iloc[0].tolist()]
        header_tokens = {
            "barcode": ("barcode", "bc", "barcode_id", "cell_id", "id"),
            "in_tissue": ("in_tissue", "tissue", "in_tissue_flag"),
            "x_pxl": ("pxl_col_in_fullres", "pxl_col", "x_pixel", "xpxl", "x"),
            "y_pxl": ("pxl_row_in_fullres", "pxl_row", "y_pixel", "ypxl", "y"),
            "x_arr": ("array_col", "col", "x_array", "array_x"),
            "y_arr": ("array_row", "row", "y_array", "array_y"),
        }
        # 若首行命中任一已知列名，或首行是数字列名(pandas to_csv 输出)，则视为表头
        has_header = any(
            any(tok in first_row for tok in tokens)
            for tokens in header_tokens.values()
        ) or all(tok.isdigit() for tok in first_row)
        if has_header:
            header = first_row
            data = raw.iloc[1:].reset_index(drop=True)
        else:
            header = [f"col_{i}" for i in range(raw.shape[1])]
            data = raw

        def find_col(tokens):
            """返回匹配 tokens 的列索引；优先精确匹配，其次子串匹配。"""
            for i, h in enumerate(header):
                if h in tokens:
                    return i
            for i, h in enumerate(header):
                for tok in tokens:
                    if tok in h:
                        return i
            return None

        col_bc = find_col(header_tokens["barcode"])
        col_in = find_col(header_tokens["in_tissue"])
        col_xp = find_col(header_tokens["x_pxl"])
        col_yp = find_col(header_tokens["y_pxl"])
        col_xa = find_col(header_tokens["x_arr"])
        col_ya = find_col(header_tokens["y_arr"])

        # 无表头且无任何已知列名时: 若全部为数值列, 用前两列作 x/y (唯一无表头 fallback)
        if col_xp is None and col_yp is None and col_xa is None and col_ya is None:
            numeric_cols = [i for i in range(data.shape[1])
                            if pd.api.types.is_numeric_dtype(data.iloc[:, i])]
            if len(numeric_cols) >= 2:
                col_xp, col_yp = numeric_cols[0], numeric_cols[1]

        # 像素坐标优先；无像素坐标时退回 array 坐标
        if col_xp is not None and col_yp is not None:
            x_col, y_col = col_xp, col_yp
            coord_kind = "pixel"
        elif col_xa is not None and col_ya is not None:
            x_col, y_col = col_xa, col_ya
            coord_kind = "array"
        else:
            x_col = y_col = None
            coord_kind = None

        if coord_kind is None:
            steps.append("Warning: no recognizable x/y coordinate columns in coordinate file")
        else:
            pos = data.copy()
            if col_in is not None:
                pos = pos[pos.iloc[:, col_in].astype(str).str.strip().isin(["1", "True", "true"])]
            coords = pos.iloc[:, [x_col, y_col]].astype(float)

            if col_bc is not None:
                # 有 barcode：按 obs_names 对齐
                pos.index = pos.iloc[:, col_bc].astype(str).values
                coords.index = pos.index
                missing = set(adata.obs_names) - set(coords.index)
                if missing:
                    steps.append(f"Warning: {len(missing)} spots missing from coordinate file")
                coords = coords.loc[adata.obs_names]
            else:
                # 无 barcode：按行顺序对齐，行数必须一致
                if len(coords) != adata.n_obs:
                    steps.append(f"Warning: coordinate file has {len(coords)} rows but adata has {adata.n_obs} spots; coords NOT attached")
                    coords = None
                else:
                    coords.index = adata.obs_names

            if coords is not None and len(coords) > 0:
                adata.obsm["spatial"] = coords.values.astype(float)
                steps.append(f"Attached {coord_kind} coordinates from {os.path.basename(coordinates_file)} ({len(coords)} spots)")

    elif "spatial" not in adata.obsm:

        steps.append("Warning: no coordinates attached (no coordinates_file given and h5ad lacks obsm['spatial'])")



    if sample_id:

        adata.obs["sample_id"] = sample_id

        steps.append(f"Set sample_id = {sample_id}")



    adata.uns["spatial_meta"] = {

        "source": counts_file,

        "n_spots": adata.n_obs,

        "n_genes": adata.n_vars,

        "has_coords": "spatial" in adata.obsm,

    }



    adata.write(output_path, compression="lzf")

    steps.append(f"Saved loaded AnnData to {output_path}")



    # Plot spot layout

    if "spatial" in adata.obsm:

        _plot_setup(os.path.dirname(output_path) or ".")

        import matplotlib.pyplot as plt



        if plot_path is None:

            plot_path = os.path.splitext(output_path)[0] + "_layout.png"

        fig, ax = plt.subplots(1, 1, figsize=(6, 6))

        _spatial_scatter(adata, None, ax, "Visium spot layout")

        fig.savefig(plot_path, dpi=150, bbox_inches="tight")

        plt.close(fig)

        steps.append(f"Spot layout plot saved to {plot_path}")

    else:

        steps.append("Warning: no coordinates to plot; spot layout plot skipped")



    return "\n".join(steps)





def filter_visium_spots(

    adata_path: str,

    output_path: str = "./visium_filtered.h5ad",

    plot_path: str = None,

    min_counts: int = 200,

    min_genes: int = 20,

    pct_mt: float = 20.0,

) -> str:

    """Filter low-quality spots in 10x Visium data.



    Applies standard QC filtering to a Visium AnnData: removes spots with too

    few counts or too few detected genes, and spots with a high fraction of

    mitochondrial reads. Adds QC metrics to ``adata.obs`` and plots the QC

    violin distributions before and after filtering.



    Args:

        adata_path: Path to the input AnnData (h5ad) file.

        output_path: Path where the filtered AnnData will be saved.

        plot_path: Path where the QC plot will be saved.

        min_counts: Minimum number of total counts per spot.

        min_genes: Minimum number of detected genes per spot.

        pct_mt: Maximum allowed percentage of mitochondrial counts.



    Returns:

        str: Step-by-step log describing the filtering result.

    """

    import scanpy as sc



    steps = []

    adata = sc.read_h5ad(adata_path)

    steps.append(f"Loaded AnnData: {adata.shape[0]} spots x {adata.shape[1]} genes")



    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)

    n_before = adata.n_obs



    # Plot QC metrics before filtering

    _plot_setup(os.path.dirname(output_path) or ".")

    import matplotlib.pyplot as plt



    if plot_path is None:

        plot_path = os.path.splitext(output_path)[0] + "_qc.png"

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    for ax_i, key in enumerate(["n_genes_by_counts", "total_counts", "pct_counts_mt"]):

        sc.pl.violin(adata, key, jitter=0.4, ax=axes[ax_i], show=False)

    fig.suptitle("QC metrics (before filtering)")

    fig.savefig(plot_path, dpi=150, bbox_inches="tight")

    plt.close(fig)

    steps.append(f"QC metrics plot saved to {plot_path}")



    adata = adata[adata.obs["total_counts"] >= min_counts, :].copy()

    adata = adata[adata.obs["n_genes_by_counts"] >= min_genes, :].copy()

    adata = adata[adata.obs["pct_counts_mt"] <= pct_mt, :].copy()



    n_after = adata.n_obs

    steps.append(

        f"QC filtering: kept {n_after}/{n_before} spots "

        f"(min_counts={min_counts}, min_genes={min_genes}, pct_mt<={pct_mt})"

    )



    adata.write(output_path, compression="lzf")

    steps.append(f"Saved filtered AnnData to {output_path}")



    return "\n".join(steps)





def normalize_visium(

    adata_path: str,

    output_path: str = "./visium_normalized.h5ad",

    plot_path: str = None,

    target_sum: float = 1e4,

    n_top_genes: int = 2000,

) -> str:

    """Normalize and log-transform 10x Visium data.



    Applies the standard scanpy pipeline: normalize_total (library size

    normalization) + log1p, optionally followed by highly variable gene

    selection. Raw counts are stored in ``adata.raw``. Plots the total-counts

    distribution before/after normalization and the HVG expression plot.



    Args:

        adata_path: Path to the input AnnData (h5ad) file.

        output_path: Path where the normalized AnnData will be saved.

        plot_path: Path where the normalization plot will be saved.

        target_sum: Target library size for normalize_total.

        n_top_genes: Number of highly variable genes to select (0 to skip).



    Returns:

        str: Step-by-step log describing the normalization result.

    """

    import scanpy as sc



    steps = []

    adata = sc.read_h5ad(adata_path)

    steps.append(f"Loaded AnnData: {adata.shape[0]} spots x {adata.shape[1]} genes")



    # Keep raw counts in .raw for downstream deconvolution tools

    adata = _raw_count_adata(adata, "normalization input")
    adata.layers["counts"] = adata.X.copy()
    adata.raw = adata.copy()

    steps.append("Stored raw counts in adata.raw")



    sc.pp.normalize_total(adata, target_sum=target_sum)

    sc.pp.log1p(adata)

    steps.append(f"Normalized to target_sum={target_sum} and log1p-transformed")



    if n_top_genes and n_top_genes > 0:

        # HVG on raw counts (seurat_v3 flavor requires raw count data)

        raw_adata = adata.raw.to_adata()

        flavor = "seurat_v3"
        try:
            sc.pp.highly_variable_genes(raw_adata, n_top_genes=n_top_genes, flavor=flavor)
            for column in ("highly_variable", "highly_variable_rank", "means", "variances", "variances_norm"):
                adata.var[column] = raw_adata.var[column]
        except (ImportError, ValueError) as exc:
            flavor = "seurat"
            steps.append(f"seurat_v3 unavailable ({exc}); using seurat")
            sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=flavor)
        adata.uns["hvg_method"] = flavor
        n_hvg = int(adata.var["highly_variable"].sum())
        steps.append(f"Selected {n_hvg} highly variable genes ({flavor}, top {n_top_genes})")

    adata.write(output_path, compression="lzf")

    steps.append(f"Saved normalized AnnData to {output_path}")



    # Plot: HVG dispersion if HVG was computed, else total-counts histogram

    _plot_setup(os.path.dirname(output_path) or ".")

    import matplotlib.pyplot as plt



    if plot_path is None:

        plot_path = os.path.splitext(output_path)[0] + "_norm.png"

    if n_top_genes and n_top_genes > 0 and "highly_variable" in adata.var:

        fig, ax = plt.subplots(1, 1, figsize=(6, 5))

        metric = "variances_norm" if adata.uns["hvg_method"] == "seurat_v3" else "dispersions_norm"
        hvg_df = adata.var[["means", metric, "highly_variable"]].copy()

        ax.scatter(hvg_df.loc[~hvg_df["highly_variable"], "means"],

                   hvg_df.loc[~hvg_df["highly_variable"], metric],

                   s=5, c="lightgrey", label="not HVG")

        ax.scatter(hvg_df.loc[hvg_df["highly_variable"], "means"],

                   hvg_df.loc[hvg_df["highly_variable"], metric],

                   s=5, c="red", label="HVG")

        ax.set_xlabel("mean expression")

        ax.set_ylabel("normalized variance" if metric == "variances_norm" else "normalized dispersion")

        ax.set_title("Highly variable genes")

        ax.legend()

        fig.savefig(plot_path, dpi=150, bbox_inches="tight")

        plt.close(fig)

        steps.append(f"HVG plot saved to {plot_path}")

    else:

        fig, ax = plt.subplots(1, 1, figsize=(6, 5))

        ax.hist(np.asarray(adata.layers["counts"].sum(axis=1)).ravel(), bins=50, color="steelblue")

        ax.set_xlabel("total_counts (raw)")

        ax.set_ylabel("spots")

        ax.set_title("Raw library size distribution")

        fig.savefig(plot_path, dpi=150, bbox_inches="tight")

        plt.close(fig)

        steps.append(f"Normalization plot saved to {plot_path}")



    return "\n".join(steps)





# ============================================================================

# Stage 2: Clustering and spatial domain detection

# ============================================================================





def cluster_spatial_data(

    adata_path: str,

    output_path: str = "./visium_clustered.h5ad",

    plot_path: str = None,

    n_pcs: int = 30,

    resolution: float = None,

) -> str:

    """Cluster 10x Visium spots using the standard PCA + UMAP + Leiden pipeline.



    Plots the UMAP colored by Leiden clusters and the spatial scatter colored

    by clusters.



    Args:

        adata_path: Path to the normalized AnnData (h5ad) file.

        output_path: Path where the clustered AnnData will be saved.

        plot_path: Path where the cluster plots will be saved (two panels).

        n_pcs: Number of principal components to use.

        resolution: Leiden clustering resolution.



    Returns:

        str: Step-by-step log describing the clustering result.

    """

    import scanpy as sc



    steps = []

    adata = sc.read_h5ad(adata_path)

    steps.append(f"Loaded AnnData: {adata.shape[0]} spots x {adata.shape[1]} genes")



    if resolution is None:
        resolution = 1.0
    n_features = int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var else adata.n_vars
    n_pcs = min(n_pcs, adata.n_obs - 1, n_features - 1)
    if n_pcs < 1:
        raise ValueError("PCA requires at least two observations and two selected genes")
    sc.pp.pca(adata, n_comps=n_pcs)

    steps.append(f"Computed PCA with {n_pcs} components")



    sc.pp.neighbors(adata, n_pcs=n_pcs)

    sc.tl.umap(adata)

    sc.tl.leiden(adata, resolution=resolution)

    n_clusters = adata.obs["leiden"].nunique()

    steps.append(f"Leiden clustering (res={resolution}) found {n_clusters} clusters")



    adata.write(output_path, compression="lzf")

    steps.append(f"Saved clustered AnnData to {output_path}")



    # Plot UMAP + spatial

    _plot_setup(os.path.dirname(output_path) or ".")

    import matplotlib.pyplot as plt



    if plot_path is None:

        plot_path = os.path.splitext(output_path)[0] + "_clusters.png"

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    sc.pl.umap(adata, color="leiden", ax=axes[0], show=False, legend_loc="on data")

    axes[0].set_title("UMAP (Leiden clusters)")

    _spatial_scatter(adata, "leiden", axes[1], "Spatial (Leiden clusters)")

    fig.savefig(plot_path, dpi=150, bbox_inches="tight")

    plt.close(fig)

    steps.append(f"Cluster plots saved to {plot_path}")



    return "\n".join(steps)



def identify_spatial_domains(
    adata_path: str,
    output_path: str = "./visium_domains.h5ad",
    plot_path: str = None,
    resolution: float = None,
    n_neighs: int = None,
    n_pcs: int = None,
    n_top_genes: int = None,
    min_domains: int = 5,
    max_domains: int = 25,
    spatial_weight: float = 0.0,
) -> str:
    """Cluster expression, optionally blending a spatial nearest-neighbor graph.

    ``spatial_weight=0`` preserves expression-only clustering; values in (0, 1]
    mix symmetrically degree-normalized expression and spatial graphs as
    (1 - weight) * expression + weight * spatial. Coordinates must share a
    common physical scale across axes. The method is recorded in ``uns``.

    FIX 2026-08-26: 原实现只用空间坐标 6-近邻图 + 固定 resolution=1.0 跑 Leiden,
    在分割后细胞级数据 (数十万细胞) 上会严重过分割 (实测 308996 细胞 -> 237 个域).
    该设计不可泛化: 不同平台 (spot/bin/cell) 的图规模差异巨大, 固定参数必然失效.

    新实现 (scanpy 标准流程, 数据驱动):
      1. 基因表达驱动 (scanpy 官方标准流程): HVG -> normalize+log1p(若 counts)
         -> scale(小数据) -> PCA -> KNN 邻域图 -> Leiden (flavor=igraph)
         (分割后细胞级数据按 scRNA 处理, 10x 官方指南)
      2. 参数决策: resolution 默认 1.0 (官方标准); 由调用方 LLM 根据数据画像
         判断后显式传参覆盖. 不做参数自动搜索 (LLM-in-the-loop 原则 2026-08-26).
      3. 兼容 spot/bin/cell 任意规模: PCA 用 arpack (可处理稀疏大矩阵),
         大数据跳过 scale 避免稠密化 (30 万细胞实测).

    Args:
        adata_path: Path to the normalized AnnData (h5ad).
        output_path: Path where the domain-annotated AnnData will be saved.
        plot_path: Path where the domain plot will be saved.
        resolution: Leiden resolution. If None, default 1.0 (LLM may override).
        n_neighs: KNN neighbors for expression graph. If None, data-driven
                  (min(15, max(5, n_obs//1000))).
        n_pcs: Number of PCA components. If None, min(30, n_obs-1, n_vars-1).
        n_top_genes: Number of HVGs for PCA. If None, min(2000, n_vars).
        min_domains, max_domains: (kept for API compat; no auto-search).

    Returns:
        str: Step-by-step log describing the spatial domain result.
    """
    import scanpy as sc
    import numpy as np

    steps = []
    adata = sc.read_h5ad(adata_path)
    steps.append(f"Loaded AnnData: {adata.shape[0]} spots x {adata.shape[1]} genes")

    if "spatial" not in adata.obsm:
        return "Error: no spatial coordinates found in adata.obsm['spatial']; run load_visium_data first"

    # FIX 2026-08-26: 数据规模自适应默认 (能算的直接算, LLM 可显式覆盖)
    n_obs, n_vars = adata.shape
    if n_obs < 3 or n_vars < 2:
        raise ValueError("Domain clustering requires at least 3 observations and 2 genes")
    if not np.isfinite(spatial_weight) or not 0 <= spatial_weight <= 1:
        raise ValueError("spatial_weight must be between 0 and 1")
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    if n_top_genes is None:
        n_top_genes = min(2000, n_vars)
    if n_pcs is None:
        n_pcs = min(30, n_obs - 1, n_vars - 1)
    if n_neighs is None:
        n_neighs = min(15, max(5, n_obs // 1000))
    n_neighs = min(n_obs - 1, max(2, n_neighs))
    steps.append(
        f"Data-driven params: n_top_genes={n_top_genes}, n_pcs={n_pcs}, "
        f"k={n_neighs} (LLM may override)")

    # 1. 基因表达驱动: HVG -> scale -> PCA
    if "spatial_domain" in adata.obs.columns:
        adata.obs = adata.obs.drop(columns=["spatial_domain"])

    # FIX 2026-08-26: 输入数据可能是 counts (整数, max 可达数百) 或 log-normalized.
    # scanpy flavor="seurat" 会对输入做 expm1 (假设 log 数据), counts 会溢出 inf 报错.
    # 泛化方案: 自动检测 counts -> 用 flavor="seurat_v3" (专为 counts 设计);
    # log 数据 -> 用 flavor="seurat". 同时防御非有限值.
    import numpy as _np
    _x = adata.X
    from scipy.sparse import issparse
    if issparse(_x):
        _d = _x.data
        _is_counts = _d.size > 0 and bool(_np.allclose(_d, _np.round(_d))) and _d.min() >= 0
    else:
        _is_counts = bool(_np.allclose(_x, _np.round(_x))) and _x.min() >= 0
    # counts -> flavor="cell_ranger" (无额外依赖, scanpy 内置均值-方差归一化);
    # log 数据 -> flavor="seurat" (依赖 scikit-misc 的 seurat_v3 会额外装包, 不用).
    _flavor = "cell_ranger" if _is_counts else "seurat"
    print(f"  HVG flavor: {_flavor} (data {'counts' if _is_counts else 'log'})")

    if adata.n_vars > 500:
        sc.pp.highly_variable_genes(adata, n_top_genes=min(n_top_genes, adata.n_vars),
                                    flavor=_flavor, subset=False)
        hv = adata.var["highly_variable"].values
        adata_tmp = adata[:, hv].copy()
        steps.append(f"Selected {int(hv.sum())} HVGs")
    else:
        adata_tmp = adata.copy()
        steps.append("Few genes, skipped HVG selection")

    # FIX 2026-08-26: sc.pp.scale 会把稀疏矩阵稠密化 (30万x2000 = 6亿元素),
    # 导致 30 万细胞数据卡死. 对大规模数据跳过 scale (counts 用 log1p 后 PCA 已足够,
    # scanpy 社区标准: normalize_total + log1p 后直接 PCA).
    # 小数据 (<10万) 仍做 scale 保证与标准流程一致.
    # 另外: 输入可能是 counts (整数), PCA 前需 log1p (scanpy 对非 log 数据的要求).
    if _is_counts:
        sc.pp.normalize_total(adata_tmp, target_sum=1e4)
        sc.pp.log1p(adata_tmp)
        steps.append("Normalized + log1p counts before PCA")
    if adata_tmp.shape[0] <= 100000:
        sc.pp.scale(adata_tmp, max_value=10)
    else:
        steps.append("Skip scale for large data (avoid densification)")
    n_comps = min(n_pcs, adata_tmp.shape[0] - 1, adata_tmp.shape[1] - 1)
    sc.tl.pca(adata_tmp, n_comps=n_comps, svd_solver="arpack")
    adata.obsm["X_pca_domain"] = adata_tmp.obsm["X_pca"]
    steps.append(f"PCA ({n_comps} components) on expression")

    # 2. 表达邻域图 (scanpy 标准 KNN, k 来自数据规模自适应)
    sc.pp.neighbors(adata, n_neighbors=n_neighs, n_pcs=n_comps,
                    use_rep="X_pca_domain", key_added="domain_neighbors")
    steps.append(f"Expression KNN graph (k={n_neighs})")

    method = "expression"
    if spatial_weight > 0:
        from scipy.sparse import diags
        from sklearn.neighbors import kneighbors_graph
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        if coords.shape != (n_obs, 2) or not np.isfinite(coords).all():
            raise ValueError("Spatial graph requires finite (n_obs, 2) coordinates")
        spatial_graph = kneighbors_graph(coords, n_neighs, mode="connectivity",
                                         include_self=False)
        spatial_graph = spatial_graph.maximum(spatial_graph.T).tocsr()
        def normalize_graph(graph):
            degree = np.asarray(graph.sum(axis=1)).ravel()
            inv = np.zeros_like(degree, dtype=float)
            np.divide(1.0, np.sqrt(degree), out=inv, where=degree > 0)
            d = diags(inv)
            return (d @ graph @ d).tocsr()
        graph_key = adata.uns["domain_neighbors"]["connectivities_key"]
        adata.obsp["domain_expression_connectivities"] = adata.obsp[graph_key].copy()
        adata.obsp["domain_spatial_connectivities"] = spatial_graph
        adata.obsp[graph_key] = (
            (1 - spatial_weight) * normalize_graph(adata.obsp[graph_key])
            + spatial_weight * normalize_graph(spatial_graph))
        method = "expression_spatial"
        steps.append(f"Blended expression/spatial graph (spatial_weight={spatial_weight})")
    adata.uns["domain_clustering"] = {"method": method,
                                      "spatial_weight": float(spatial_weight),
                                      "n_neighbors": int(n_neighs)}

    # 3. Leiden (官方标准: 默认 resolution=1.0; 用户显式传参才覆盖)
    if resolution is None:
        resolution = 1.0
    # FIX 2026-08-26: 默认 leidenalg 后端在 30 万细胞图上卡死 (30+ min 无结果),
    # igraph flavor 25 秒完成 (scanpy 官方推荐的新后端, 见 FutureWarning).
    sc.tl.leiden(adata, resolution=resolution, key_added="spatial_domain",
                 neighbors_key="domain_neighbors", flavor="igraph",
                 n_iterations=2, directed=False)
    n_domains = adata.obs["spatial_domain"].nunique()
    steps.append(f"{method} clustering found {n_domains} groups (res={resolution})")

    # 5. Domain 注释 (聚类后必须给每个 domain 生物学证据, 供 LLM 按命名规范命名)
    #    代码产出证据 (annotation CSV): marker 基因 + 细胞类型组成.
    #    LLM 命名流程 (5 步, 见 spatial-visium-hd/SKILL.md "空间域命名规范"):
    #      1. 读 annotation CSV (marker 基因 + 细胞类型组成)
    #      2. 判定组织类型: 从样本元数据/组织来源推断 (肠/肺/心/肝/癌)
    #      3. 对照组织类型参考表选候选名: 优先解剖结构, 其次组织学特征
    #      4. 用 marker 基因验证: 候选名必须有 marker 支持
    #         (如"隐窝底部"需 Lgr5/Defa 等潘氏/干细胞 marker; "心肌层"需 TNNT2/ACTA2 等)
    #      5. 输出统一格式: D{n} {结构名}({特征}); 证据不足 → D{n} 未命名(待定) 并说明缺什么证据
    #    命名不在代码里写死 (避免 LLM 自由发挥, 统一走 5 步流程).
    try:
        if adata.n_vars <= 100000:
            # 30 万细胞 rank_genes_groups 很慢, 用 downsample (每 domain 最多 5000)
            n_per_dom = adata.obs.groupby("spatial_domain", observed=True).size()
            max_per = int(n_per_dom.min())
            sub_idx = adata.obs.groupby("spatial_domain", observed=True).apply(
                lambda x: x.sample(min(len(x), 5000), random_state=0)
            ).index.get_level_values(-1)
            adata_sub = adata[sub_idx].copy()
            if _is_counts:
                sc.pp.normalize_total(adata_sub, target_sum=1e4)
                sc.pp.log1p(adata_sub)
            sc.tl.rank_genes_groups(adata_sub, "spatial_domain", method="wilcoxon",
                                    n_genes=10, use_raw=False)
            # 汇总: 每 domain top marker + 细胞类型组成
            import pandas as pd
            rows = []
            dom_cats = list(adata.obs["spatial_domain"].astype("category").cat.categories)
            for d in dom_cats:
                # marker genes
                try:
                    de = adata_sub.uns["rank_genes_groups"]
                    idx = list(de["names"].dtype.names).index(str(d))
                    genes = list(de["names"][str(d)][:10])
                    lfc = [round(float(x), 3) for x in de["logfoldchanges"][str(d)][:10]]
                    marker_str = "; ".join(f"{g}({l})" for g, l in zip(genes, lfc))
                except Exception:
                    marker_str = ""
                # cell type composition
                if "cell_type" in adata.obs.columns:
                    comp = adata.obs.loc[adata.obs["spatial_domain"] == d, "cell_type"] \
                        .value_counts(normalize=True).head(3)
                    comp_str = "; ".join(f"{k}:{v:.0%}" for k, v in comp.items())
                else:
                    comp_str = "N/A (no cell_type)"
                rows.append({"domain": d, "n_cells": int(n_per_dom.get(d, 0)),
                             "top_markers": marker_str, "cell_types": comp_str})
            df = pd.DataFrame(rows)
            ann_path = os.path.splitext(output_path)[0] + "_annotation.csv"
            df.to_csv(ann_path, index=False)
            steps.append(f"Domain annotation saved to {ann_path}")
            steps.append("Domain annotation (LLM uses these to name domains):")
            for _, r in df.iterrows():
                steps.append(f"  D{r['domain']}: n={r['n_cells']} | "
                             f"markers={r['top_markers'][:80]} | "
                             f"cells={r['cell_types'][:60]}")
        else:
            steps.append("Skip domain annotation (>100k vars, too heavy)")
    except Exception as e:
        steps.append(f"Domain annotation skipped (error: {str(e)[:100]})")

    adata.write(output_path, compression="lzf")
    steps.append(f"Saved domain-annotated AnnData to {output_path}")

    # 4. Plot domains on spatial coordinates
    _plot_setup(os.path.dirname(output_path) or ".")
    import matplotlib.pyplot as plt

    if plot_path is None:
        plot_path = os.path.splitext(output_path)[0] + "_domains.png"
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    _spatial_scatter(adata, "spatial_domain", ax, f"Domains ({method})")
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    steps.append(f"Spatial domain plot saved to {plot_path}")

    return "\n".join(steps)


def find_spatially_variable_genes(

    adata_path: str,

    output_path: str = "./svg_results.csv",

    plot_path: str = None,

    n_top_genes: int = None,

    genes_to_plot: list = None,

    n_jobs: int = 1,

) -> str:

    """Detect spatially variable genes (SVGs) in 10x Visium data using Moran's I.



    Uses squidpy's spatial autocorrelation (Moran's I) to rank genes by spatial

    clustering of expression. Saves the ranked SVG table to CSV and plots the

    top SVG expressions on the spatial coordinates.



    Args:

        adata_path: Path to the normalized AnnData (h5ad) file with spatial coordinates.

        output_path: Path where the SVG ranking CSV will be saved.

        plot_path: Path where the top-SVG spatial plots will be saved.

        n_top_genes: Number of top spatially variable genes to report.

        genes_to_plot: Optional list of gene names to plot on spatial coordinates.

            If provided, these genes are plotted (must exist in var_names).

            If None, the top 4 SVGs are plotted by default.

        n_jobs: Number of parallel jobs for autocorrelation computation.



    Returns:

        str: Step-by-step log describing the SVG detection result.

    """

    try:

        import squidpy as sq

    except ImportError:

        return "Error: squidpy is required for SVG detection. Install with: pip install squidpy"



    steps = []

    adata = sc.read_h5ad(adata_path)

    steps.append(f"Loaded AnnData: {adata.shape[0]} spots x {adata.shape[1]} genes")



    if "spatial" not in adata.obsm:

        return "Error: no spatial coordinates found in adata.obsm['spatial']; run load_visium_data first"

    # FIX 2026-08-26: n_top_genes 数据驱动 (LLM 可覆盖)
    if n_top_genes is None:
        n_top_genes = min(100, adata.n_vars)

    if "spatial_connectivities" not in adata.obsp:
        _k = min(15, max(5, adata.shape[0] // 1000))
        sq.gr.spatial_neighbors(adata, n_neighs=_k, coord_type="generic")
        steps.append(f"Built spatial neighborhood graph (k={_k})")

    _n_genes = min(n_top_genes, adata.n_vars)
    sq.gr.spatial_autocorr(

        adata,

        mode="moran",

        n_jobs=n_jobs,

        genes=adata.var_names,

    )

    steps.append("Computed Moran's I for all genes")



    moran = adata.uns["moranI"]

    moran = moran.sort_values("I", ascending=False)

    top = moran.head(n_top_genes)

    top.to_csv(output_path)

    steps.append(f"Top {len(top)} spatially variable genes saved to {output_path}")



    # Plot top 4 SVGs on spatial coordinates

    _plot_setup(os.path.dirname(output_path) or ".")

    import matplotlib.pyplot as plt



    if plot_path is None:

        plot_path = os.path.splitext(output_path)[0] + "_top_svgs.png"

    if genes_to_plot:

        missing = [g for g in genes_to_plot if g not in adata.var_names]

        if missing:

            return "\\n".join(steps) + f"\\nError: genes not found in var_names: {missing}"

        plot_genes = genes_to_plot

        title = "Spatially variable genes (user selected)"

    else:

        plot_genes = top.index[:4].tolist()

        title = "Top spatially variable genes"

    n_plots = len(plot_genes)

    ncols = min(n_plots, 4)

    nrows = (n_plots + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows), squeeze=False)

    for ax, gene in zip(axes.ravel(), plot_genes):

        moran_i = top.loc[gene, "I"] if gene in top.index else float("nan")

        _spatial_scatter(adata, gene, ax, f"{gene} (Moran's I={moran_i:.3f})" if gene in top.index else gene)

    for ax in axes.ravel()[len(plot_genes):]:

        ax.axis("off")

    fig.suptitle(title)

    fig.tight_layout()

    fig.savefig(plot_path, dpi=150, bbox_inches="tight")

    plt.close(fig)

    steps.append(f"SVG spatial plots ({len(plot_genes)} genes) saved to {plot_path}")



    return "\n".join(steps)





# ============================================================================

# Stage 4: Cell-type deconvolution

# ============================================================================





def _plot_proportions(prop_df, plot_path: str, title: str) -> None:

    """Plot per-spot cell-type proportions as a stacked bar/area summary."""

    _plot_setup(os.path.dirname(plot_path) or ".")

    import matplotlib.pyplot as plt



    # Bar chart of mean proportions per cell type (numeric columns only)

    num_cols = prop_df.select_dtypes(include="number").columns

    mean_prop = prop_df[num_cols].mean(axis=0).sort_values(ascending=False)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].bar(range(len(mean_prop)), mean_prop.values, color="steelblue")

    axes[0].set_xticks(range(len(mean_prop)))

    axes[0].set_xticklabels(mean_prop.index, rotation=90, fontsize=7)

    axes[0].set_ylabel("mean proportion")

    axes[0].set_title("Mean cell-type proportions")



    # Stacked proportions for first 20 spots (sampled)

    # Each bar = one spot, stack = cell-type proportions (sums to 1)

    n_show = min(20, prop_df.shape[0])

    sample = prop_df[num_cols].iloc[:n_show]

    sample.plot(kind="bar", stacked=True, ax=axes[1], colormap="tab20", width=0.9)

    axes[1].set_xlabel("spot")

    axes[1].set_ylabel("proportion")

    axes[1].set_title(f"Per-spot proportions (first {n_show} spots)")

    axes[1].legend(fontsize=5, ncol=2)

    fig.suptitle(title)

    fig.tight_layout()

    fig.savefig(plot_path, dpi=150, bbox_inches="tight")

    plt.close(fig)





def deconvolve_spatial_spotlight(

    ref_h5ad: str,

    st_h5ad: str,

    cell_type_key: str,

    output_dir: str = "./spotlight_results",

    plot_path: str = None,

    r_script_path: str = "Rscript",

) -> str:

    """Deconvolve 10x Visium spots into cell-type proportions using SPOTlight.



    SPOTlight (NAR 2021) performs seeded NMF on reference scRNA-seq data to

    learn cell-type topic profiles, then deconvolves each spatial spot into

    cell-type proportions via non-negative least squares. This tool bridges to

    the R implementation through a temporary R script. Saves a proportions

    summary plot.



    Args:

        ref_h5ad: Path to the reference scRNA-seq AnnData (h5ad) with cell-type labels.

        st_h5ad: Path to the Visium AnnData (h5ad) to deconvolve.

        cell_type_key: Column in ref_h5ad.obs containing cell-type labels.

        output_dir: Directory where SPOTlight results (proportions CSV) are saved.

        plot_path: Path where the proportions plot will be saved.

        r_script_path: Path to the Rscript executable.



    Returns:

        str: Step-by-step log describing the deconvolution result.

    """

    import subprocess

    import tempfile



    if not os.path.exists(ref_h5ad) or not os.path.exists(st_h5ad):

        return "Error: ref_h5ad or st_h5ad file not found"



    os.makedirs(output_dir, exist_ok=True)

    steps = []

    steps.append(f"SPOTlight deconvolution: ref={ref_h5ad}, st={st_h5ad}")



    # Convert h5ad inputs to 10x format for R (deterministic, no R h5ad reader needed)

    import scipy.io as sio



    bridge_dir = os.path.join(output_dir, "_bridge")

    os.makedirs(bridge_dir, exist_ok=True)



    def _h5ad_to_10x(h5ad_path, prefix):

        """Write an h5ad to 10x mtx/features/barcodes files (genes x cells).



        Uses gene symbols when available (var['feature_name']), falling back

        to var_names (Ensembl IDs). SPOTlight/scran gene filters need symbols.

        """

        ad = _raw_count_adata(sc.read_h5ad(h5ad_path), "SPOTlight input")

        # Gene symbols for R-side filters (^Rp[l|s]|Mt etc.)

        if "feature_name" in ad.var.columns:

            symbols = ad.var["feature_name"].astype(str).values

            # Fill missing/empty symbols with var_names

            symbols = np.where((symbols == "") | (symbols == "nan"), ad.var_names.astype(str), symbols)

        else:

            symbols = ad.var_names.astype(str).values

        mtx_path = os.path.join(bridge_dir, f"{prefix}_matrix.mtx")

        feat_path = os.path.join(bridge_dir, f"{prefix}_features.tsv")

        bar_path = os.path.join(bridge_dir, f"{prefix}_barcodes.tsv")

        from scipy.sparse import csr_matrix
        sio.mmwrite(mtx_path, csr_matrix(ad.X.T))

        with open(feat_path, "w") as f:

            for g in symbols:

                f.write(f"{g}\t{g}\tGene Expression\n")

        with open(bar_path, "w") as f:

            for b in ad.obs_names:

                f.write(f"{b}\n")

        if cell_type_key in ad.obs:

            ct_path = os.path.join(bridge_dir, f"{prefix}_celltypes.tsv")

            ad.obs[[cell_type_key]].to_csv(ct_path, header=True, index=False)

        return mtx_path, feat_path, bar_path



    steps.append("Converting h5ad to 10x format for R bridge...")

    ref_mtx, ref_feat, ref_bar = _h5ad_to_10x(ref_h5ad, "ref")

    st_mtx, st_feat, st_bar = _h5ad_to_10x(st_h5ad, "st")

    steps.append("Conversion done")



    # Write the R bridge script

    r_script = f"""

suppressPackageStartupMessages({{

  library(SPOTlight)

  library(SingleCellExperiment)

  library(scater)

  library(scran)

  library(Matrix)

}})



read10x <- function(mtx, feat, bar) {{

  cnt <- readMM(mtx)  # genes x cells (10x standard)

  rownames(cnt) <- read.delim(feat, header = FALSE)[,1]   # genes

  colnames(cnt) <- read.delim(bar, header = FALSE)[,1]    # cells

  sce <- SingleCellExperiment(assays = list(counts = cnt))

  return(sce)

}}



ref <- read10x({_r_literal(ref_mtx)}, {_r_literal(ref_feat)}, {_r_literal(ref_bar)})

st <- read10x({_r_literal(st_mtx)}, {_r_literal(st_feat)}, {_r_literal(st_bar)})



# Cell-type labels

ct <- read.delim({_r_literal(os.path.join(bridge_dir, "ref_celltypes.tsv"))})

colData(ref)$cell_type <- ct[,1]



# Downsample to max 100 cells per type (official recommendation)

set.seed(123)

idx <- split(seq(ncol(ref)), ref$cell_type)

cs_keep <- lapply(idx, function(i) {{

  sample(i, min(length(i), 100))

}})

ref <- ref[, unlist(cs_keep)]

ref <- logNormCounts(ref)



# HVGs

genes <- !grepl("^Rp[l|s]|Mt", rownames(ref))

dec <- modelGeneVar(ref, subset.row = genes)

hvg <- getTopHVGs(dec, n = 3000)

colLabels(ref) <- colData(ref)$cell_type



# Marker genes -> data.frame format (SPOTlight 1.10 official tutorial)

# Keep top n_top markers per cell type (default 50) to keep NMF init fast.

mgs <- scoreMarkers(ref, subset.row = genes)

mgs_fil <- lapply(names(mgs), function(i) {{

  x <- mgs[[i]]

  x <- x[order(x$mean.AUC, decreasing = TRUE), ]

  x <- head(x, 50)

  x$gene <- rownames(x)

  x$cluster <- i

  data.frame(x)

}})

mgs_df <- do.call(rbind, mgs_fil)

cat("mgs_df:", nrow(mgs_df), "rows\n")



# SPOTlight (v1.10: x=SCE, y=spatial SCE, mgs=data.frame)

res <- SPOTlight(

  x = ref,

  y = st,

  groups = as.character(ref$cell_type),

  mgs = mgs_df,

  hvg = hvg,

  weight_id = "mean.AUC",

  group_id = "cluster",

  gene_id = "gene"

)

prop <- as.data.frame(res$mat)

prop$barcode <- rownames(prop)

write.csv(prop, file.path({_r_literal(output_dir)}, 'spotlight_proportions.csv'), row.names = FALSE)

cat("SPOTlight done\\n")

"""

    with tempfile.NamedTemporaryFile(suffix=".R", mode="w", encoding="utf-8", delete=False, dir=output_dir) as f:

        f.write(r_script)

        r_script_path_tmp = f.name



    try:

        result = _run_r_script(r_script_path, r_script_path_tmp, output_dir)

        steps.append(f"R output: {result.stdout.strip()[:500]}")

        if result.returncode != 0:

            steps.append(f"R error: {result.stderr.strip()[:500]}")

            return "\n".join(steps)

        out_csv = os.path.join(output_dir, "spotlight_proportions.csv")

        if os.path.exists(out_csv):

            prop = pd.read_csv(out_csv).set_index("barcode")

            steps.append(f"SPOTlight proportions: {prop.shape[0]} spots x {prop.shape[1]} cell types")

            steps.append(f"Results saved to {out_csv}")

            if plot_path is None:

                plot_path = os.path.join(output_dir, "spotlight_proportions.png")

            _plot_proportions(prop, plot_path, "SPOTlight deconvolution")

            steps.append(f"Proportions plot saved to {plot_path}")

        else:

            steps.append("Warning: output CSV not found")

    except Exception as e:

        steps.append(f"Error running SPOTlight: {e}")

    finally:

        os.unlink(r_script_path_tmp)



    return "\n".join(steps)





def deconvolve_spatial_destvi(

    ref_h5ad: str,

    st_h5ad: str,

    cell_type_key: str,

    output_dir: str = "./destvi_results",

    plot_path: str = None,

    n_latent: int = 30,

    max_epochs: int = 400,

    destvi_max_epochs: int = 500,

) -> str:

    """Deconvolve 10x Visium spots into cell-type proportions using DestVI.



    DestVI (Lopez et al., Nat Biotechnol 2022) is a deep generative model from

    the scvi-tools ecosystem. It first trains a CondSCVI model on the reference

    scRNA-seq data, then trains a spatial model on the Visium data to estimate

    cell-type proportions per spot. Saves a proportions summary plot.



    Args:

        ref_h5ad: Path to the reference scRNA-seq AnnData (h5ad) with cell-type labels.

        st_h5ad: Path to the Visium AnnData (h5ad) to deconvolve.

        cell_type_key: Column in ref_h5ad.obs containing cell-type labels.

        output_dir: Directory where DestVI results (proportions CSV) are saved.

        plot_path: Path where the proportions plot will be saved.

        n_latent: Latent dimensionality for CondSCVI.

        max_epochs: Number of training epochs for the CondSCVI reference model.

        destvi_max_epochs: Number of training epochs for the DestVI spatial model.



    Returns:

        str: Step-by-step log describing the deconvolution result.

    """

    try:

        import torch

        from scvi.model import CondSCVI, DestVI

    except ImportError:

        return "Error: scvi-tools required for DestVI deconvolution. Install with: pip install scvi-tools"



    if not os.path.exists(ref_h5ad) or not os.path.exists(st_h5ad):

        return "Error: ref_h5ad or st_h5ad file not found"



    os.makedirs(output_dir, exist_ok=True)

    steps = []

    steps.append(f"DestVI deconvolution: ref={ref_h5ad}, st={st_h5ad}")



    use_gpu = torch.cuda.is_available()

    steps.append(f"CUDA available: {use_gpu}")



    # Load reference

    ref = sc.read_h5ad(ref_h5ad)

    steps.append(f"Reference: {ref.shape[0]} cells x {ref.shape[1]} genes")

    if cell_type_key not in ref.obs:

        return f"Error: cell_type_key {_r_literal(cell_type_key)} not found in reference obs"



    ref = _raw_count_adata(ref, "reference")



    # Load spatial data first so we can align genes before training CondSCVI

    st = _raw_count_adata(sc.read_h5ad(st_h5ad), "spatial")

    steps.append(f"Spatial: {st.shape[0]} spots x {st.shape[1]} genes")



    # Align genes between reference (Ensembl) and spatial (symbol)

    # Reference var_names may be Ensembl IDs with a 'feature_name' column holding symbols

    if "feature_name" in ref.var.columns:

        sym2ens = dict(zip(ref.var["feature_name"].astype(str), ref.var_names))

        st_symbols = np.array([str(g) for g in st.var_names])

        st_ens = np.array([sym2ens.get(s, "") for s in st_symbols])

        keep = np.array([e in set(ref.var_names) for e in st_ens])

        n_mapped = int(keep.sum())

        steps.append(f"Gene symbol->Ensembl mapped: {n_mapped}/{len(st_symbols)}")

        if n_mapped == 0:

            return "Error: no genes could be mapped between spatial (symbol) and reference (Ensembl)"

        st = st[:, keep].copy()

        st.var_names = st_ens[keep]

        common_genes = [g for g in ref.var_names if g in set(st.var_names)]

    else:

        common_genes = list(ref.var_names.intersection(st.var_names))

    if len(common_genes) == 0:

        return "Error: no common genes between reference and spatial data"



    # Align BOTH datasets to the same gene set (required by DestVI.from_rna_model)

    ref = ref[:, common_genes].copy()

    st = st[:, common_genes].copy()

    steps.append(f"Aligned on {len(common_genes)} common genes (ref and st)")



    CondSCVI.setup_anndata(ref, labels_key=cell_type_key)

    rna_model = CondSCVI(ref, n_latent=n_latent, prior="mog")

    steps.append(f"Training CondSCVI reference model ({max_epochs} epochs)...")

    rna_model.train(

        max_epochs=max_epochs,

        accelerator="cuda" if use_gpu else "cpu",

        devices=[0] if use_gpu else 1,

        batch_size=512,

    )

    steps.append("CondSCVI reference model trained")



    # Spatial model

    st.layers["counts"] = st.X.copy()

    DestVI.setup_anndata(st, layer="counts")

    spatial_model = DestVI.from_rna_model(st, rna_model)

    steps.append("Training DestVI spatial model...")

    spatial_model.train(max_epochs=destvi_max_epochs, accelerator="cuda" if use_gpu else "cpu", devices=[0] if use_gpu else 1)

    steps.append("DestVI spatial model trained")



    # Extract proportions (returns pd.DataFrame with cell-type columns)

    prop_df = spatial_model.get_proportions()

    prop_df = prop_df.reset_index(drop=True)

    prop_df.index = st.obs_names

    out_csv = os.path.join(output_dir, "destvi_proportions.csv")

    prop_df.to_csv(out_csv)

    steps.append(f"DestVI proportions: {prop_df.shape[0]} spots x {prop_df.shape[1]} cell types")

    steps.append(f"Results saved to {out_csv}")



    if plot_path is None:

        plot_path = os.path.join(output_dir, "destvi_proportions.png")

    _plot_proportions(prop_df, plot_path, "DestVI deconvolution")

    steps.append(f"Proportions plot saved to {plot_path}")



    return "\n".join(steps)





# ============================================================================

# Stage 5: Cell-cell communication inference (CellChat v2)

# ============================================================================





def infer_spatial_cell_communication(
    st_h5ad: str,
    cell_type_key: str,
    output_dir: str = "./cellchat_results",
    plot_path: str = None,
    r_script_path: str = "Rscript",
    species: str = None,
    mode: str = "expression",
    microns_per_coordinate: float = None,
    spatial_tolerance: float = None,
    interaction_range: float = 250.0,
    scale_distance: float = 0.01,
    contact_range: float = None,
    k_min: int = 10,
) -> str:
    """Infer CellChat v2.1.2 expression or distance-constrained communication.

    Expression mode preserves the historical computation and ignores coordinates.
    Spatial mode requires finite x/y coordinates, microns per coordinate unit,
    and a nonnegative spatial_tolerance in microns. interaction_range and optional
    contact_range are in microns. scale_distance rescales group distances for
    CellChat; choose it so the smallest scaled distance is >= 1 (CellChat checks
    this). k_min is the minimum number of nearby cell pairs per group pair.
    Contact-specific filtering is enabled when contact_range is provided.
    Both modes export CSV tables, the fitted RDS object, logs and run_settings.json.
    Raw integer counts are required in layers['counts'], raw.X, or X.
    """
    import json
    import tempfile
    import numpy as np
    import scipy.io as sio
    from scipy.sparse import csr_matrix

    if mode not in {"expression", "spatial"}:
        raise ValueError("mode must be expression or spatial")
    ad = _raw_count_adata(sc.read_h5ad(st_h5ad), "CellChat input")
    if cell_type_key not in ad.obs or ad.obs[cell_type_key].isna().any():
        raise ValueError("CellChat requires a nonmissing cell-type label for every observation")
    species = detect_species(ad) if species is None else species
    if species not in {"human", "mouse"}:
        raise ValueError("Set species to human or mouse")
    spatial_args = ""
    creation = f"cellchat <- createCellChat(object = obj, group.by = {_r_literal(cell_type_key)})"
    if mode == "spatial":
        positive = {"microns_per_coordinate": microns_per_coordinate,
                    "interaction_range": interaction_range, "scale_distance": scale_distance}
        if contact_range is not None:
            positive["contact_range"] = contact_range
        for name, value in positive.items():
            if value is None or not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive for spatial mode")
        if spatial_tolerance is None or not np.isfinite(spatial_tolerance) or spatial_tolerance < 0:
            raise ValueError("spatial_tolerance must be finite and nonnegative in microns")
        if isinstance(k_min, bool) or int(k_min) != k_min or k_min < 1:
            raise ValueError("k_min must be a positive integer")
        coords = np.asarray(ad.obsm.get("spatial", []), dtype=float)
        if coords.shape != (ad.n_obs, 2) or not np.isfinite(coords).all():
            raise ValueError("Spatial CellChat requires finite (n_cells, 2) spatial coordinates")
        spatial_args = (f", distance.use = TRUE, interaction.range = {float(interaction_range)}, "
                        f"scale.distance = {float(scale_distance)}, k.min = {int(k_min)}, "
                        f"contact.dependent = {'TRUE' if contact_range is not None else 'FALSE'}")
        if contact_range is not None:
            spatial_args += f", contact.range = {float(contact_range)}"

    os.makedirs(output_dir, exist_ok=True)
    plot_path = plot_path or os.path.join(output_dir, "cellchat_network.png")
    os.makedirs(os.path.dirname(plot_path) or ".", exist_ok=True)
    bridge_dir = os.path.join(output_dir, "_bridge")
    os.makedirs(bridge_dir, exist_ok=True)
    mtx_path = os.path.join(bridge_dir, "st_matrix.mtx")
    feat_path = os.path.join(bridge_dir, "st_features.tsv")
    bar_path = os.path.join(bridge_dir, "st_barcodes.tsv")
    ct_path = os.path.join(bridge_dir, "st_celltypes.tsv")
    coord_path = os.path.join(bridge_dir, "st_coords.csv")
    sio.mmwrite(mtx_path, csr_matrix(ad.X).T)
    pd.Series(ad.var_names).to_csv(feat_path, sep="\t", index=False, header=False)
    pd.Series(ad.obs_names).to_csv(bar_path, sep="\t", index=False, header=False)
    labels = ad.obs[cell_type_key].astype(str)
    if labels.str.strip().eq("").any():
        raise ValueError("CellChat labels must be nonempty")
    if (labels == "0").any():
        # Prefix all labels to preserve identity even for mixed labels such as 0/A.
        labels = "group_" + labels
    labels.to_frame(cell_type_key).to_csv(ct_path, index=False)
    if mode == "spatial":
        pd.DataFrame(coords, index=ad.obs_names, columns=["x", "y"]).to_csv(coord_path)
        creation = f"""
coords <- as.matrix(read.csv({_r_literal(coord_path)}, row.names = 1, check.names = FALSE))
coords <- coords[colnames(obj), , drop = FALSE]
cellchat <- createCellChat(object = obj, group.by = {_r_literal(cell_type_key)},
                          datatype = "spatial", coordinates = coords,
                          spatial.factors = data.frame(ratio = {float(microns_per_coordinate)},
                                                       tol = {float(spatial_tolerance)}))
"""
    settings = {"mode": mode, "species": species, "cell_type_key": cell_type_key,
                "microns_per_coordinate": microns_per_coordinate, "spatial_tolerance": spatial_tolerance,
                "interaction_range": interaction_range, "scale_distance": scale_distance,
                "contact_range": contact_range, "k_min": k_min,
                "label_mapping": dict(zip(ad.obs[cell_type_key].astype(str), labels))}
    with open(os.path.join(output_dir, "run_settings.json"), "w", encoding="utf-8") as handle:
        json.dump(settings, handle, ensure_ascii=False, indent=2)
    r_script = f"""
suppressPackageStartupMessages({{library(CellChat); library(Seurat); library(Matrix)}})
cnt <- readMM({_r_literal(mtx_path)})
rownames(cnt) <- read.delim({_r_literal(feat_path)}, header = FALSE, quote = "", comment.char = "")[,1]
colnames(cnt) <- read.delim({_r_literal(bar_path)}, header = FALSE, quote = "", comment.char = "")[,1]
obj <- CreateSeuratObject(counts = cnt, assay = "RNA")
obj <- NormalizeData(obj)
ct <- read.csv({_r_literal(ct_path)}, check.names = FALSE)
obj <- AddMetaData(obj, metadata = ct[,1], col.name = {_r_literal(cell_type_key)})
{creation}
cellchat@DB <- CellChatDB.{species}
cellchat <- subsetData(cellchat)
future::plan("sequential")
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1{spatial_args})
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file.path({_r_literal(output_dir)}, "cellchat.rds"))
if (sum(cellchat@net$count) > 0) {{
  df <- subsetCommunication(cellchat)
}} else {{
  df <- data.frame(source = character(), target = character(), ligand = character(),
                   receptor = character(), prob = numeric(), pval = numeric(),
                   interaction_name = character(), pathway_name = character())
}}
write.csv(df, file.path({_r_literal(output_dir)}, "cellchat_communication.csv"), row.names = FALSE)
if (length(cellchat@netP$pathways) > 0) {{
  df_path <- subsetCommunication(cellchat, slot.name = "netP")
}} else {{
  df_path <- data.frame(source = character(), target = character(), pathway_name = character(),
                        prob = numeric(), pval = numeric())
}}
write.csv(df_path, file.path({_r_literal(output_dir)}, "cellchat_pathways.csv"), row.names = FALSE)
png({_r_literal(plot_path)}, width = 10, height = 10, units = "in", res = 150)
if (sum(cellchat@net$count) > 0) {{
  netVisual_circle(cellchat@net$count, vertex.weight = table(cellchat@idents),
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = "CellChat ({mode})", vertex.label.cex = 1.2)
}} else {{
  plot.new(); title("CellChat ({mode}): no significant interactions")
}}
dev.off()
cat("CellChat done\\n")
"""
    with tempfile.NamedTemporaryFile(suffix=".R", mode="w", encoding="utf-8",
                                     delete=False, dir=output_dir) as handle:
        handle.write(r_script)
        script_path = handle.name
    try:
        result = _run_r_script(r_script_path, script_path, output_dir)
        if result.returncode:
            raise RuntimeError(f"CellChat failed; see {output_dir}/r_stderr.log: {result.stderr[-1500:]}")
    finally:
        os.unlink(script_path)
    comm_csv = os.path.join(output_dir, "cellchat_communication.csv")
    comm = pd.read_csv(comm_csv)
    return (f"CellChat mode: {mode}; species: {species}\n"
            f"CellChat communication pairs: {len(comm)}\nResults saved to {comm_csv}\n"
            f"Network plot saved to {plot_path}\nModel saved to {output_dir}/cellchat.rds")


def spatial_neighborhood_enrichment(

    adata_path: str,

    cell_type_key: str,

    output_path: str = "./neighborhood_enrichment.csv",

    plot_path: str = None,

    n_neighs: int = None,

) -> str:

    """Compute cell-type neighborhood enrichment in 10x Visium data.



    Uses squidpy's neighborhood enrichment analysis to test whether pairs of

    cell types (or spatial domains) are found together more often than expected

    by chance in the spatial neighborhood graph. Saves the enrichment table and

    plots the enrichment z-score heatmap.



    Args:

        adata_path: Path to the AnnData (h5ad) with spatial coordinates and

            cell-type labels in obs.

        cell_type_key: Column in adata.obs containing cell-type labels.

        output_path: Path where the enrichment CSV will be saved.

        plot_path: Path where the enrichment heatmap will be saved.

        n_neighs: Number of spatial neighbors per spot.



    Returns:

        str: Step-by-step log describing the neighborhood enrichment result.

    """

    try:

        import squidpy as sq

    except ImportError:

        return "Error: squidpy is required for neighborhood enrichment. Install with: pip install squidpy"



    steps = []

    adata = sc.read_h5ad(adata_path)

    steps.append(f"Loaded AnnData: {adata.shape[0]} spots x {adata.shape[1]} genes")



    if cell_type_key not in adata.obs:

        return f"Error: cell_type_key {_r_literal(cell_type_key)} not found in obs"

    if "spatial" not in adata.obsm:

        return "Error: no spatial coordinates found; run load_visium_data first"

    # FIX 2026-08-26: n_neighs 数据驱动 (LLM 可覆盖)
    if n_neighs is None:
        n_neighs = min(15, max(5, adata.shape[0] // 1000))

    sq.gr.spatial_neighbors(adata, n_neighs=n_neighs, coord_type="generic")

    steps.append(f"Built spatial neighborhood graph with {n_neighs} neighbors per spot")



    sq.gr.nhood_enrichment(adata, cluster_key=cell_type_key)

    steps.append("Computed neighborhood enrichment")



    enrich = adata.uns[f"{cell_type_key}_nhood_enrichment"]

    cats = adata.obs[cell_type_key].cat.categories

    df = pd.DataFrame(enrich["zscore"], index=cats, columns=cats)

    df.to_csv(output_path)

    steps.append(f"Neighborhood enrichment z-scores saved to {output_path}")



    # Plot heatmap

    _plot_setup(os.path.dirname(output_path) or ".")

    import matplotlib.pyplot as plt



    if plot_path is None:

        plot_path = os.path.splitext(output_path)[0] + "_heatmap.png"

    fig, ax = plt.subplots(1, 1, figsize=(max(6, len(cats) * 0.6), max(5, len(cats) * 0.5)))

    im = ax.imshow(df.values, cmap="RdBu_r", vmin=-max(abs(df.values.max()), abs(df.values.min())),

                   vmax=max(abs(df.values.max()), abs(df.values.min())))

    # Readable labels: numeric cluster ids become "Cluster N"; real cell-type
    # names (e.g. "T cell", "Fibroblast") are kept as-is so the heatmap carries
    # biological meaning.
    def _readable(c):

        s = str(c)

        return f"Cluster {s}" if s.isdigit() else s

    labels = [_readable(c) for c in cats]

    ax.set_xticks(range(len(cats)))

    ax.set_yticks(range(len(cats)))

    ax.set_xticklabels(labels, rotation=90, fontsize=7)

    ax.set_yticklabels(labels, fontsize=7)

    ax.set_xlabel("Cell type B")

    ax.set_ylabel("Cell type A")

    plt.colorbar(im, ax=ax, fraction=0.046, label="z-score")

    ax.set_title("Neighborhood enrichment (z-score)")

    fig.tight_layout()

    fig.savefig(plot_path, dpi=150, bbox_inches="tight")

    plt.close(fig)

    steps.append(f"Enrichment heatmap saved to {plot_path}")



    return "\n".join(steps)





# ============================================================================

# Project scaffolding

# ============================================================================





def init_spatial_project(

    project_dir: str,

    overwrite: bool = False,

) -> str:

    """Create a standardized spatial transcriptomics project directory layout.



    All spatial tools in this module expect results to live under a project

    directory with a consistent structure, so multi-sample / multi-step

    analyses stay organized:



        <project_dir>/

        âââ raw_data/

        â?  âââ visium/          # raw Visium input (h5ad, 10x dirs, CSVs)

        â?  âââ scRNA/           # reference scRNA-seq (h5ad) for deconvolution

        âââ results/

            âââ 01_loading/      # standardized h5ad after load_visium_data

            âââ 02_qc/           # QC-filtered h5ad + plots

            âââ 03_normalization/

            âââ 04_clustering/   # clusters / spatial domains

            âââ 05_deconvolution/# SPOTlight / DestVI proportions

            âââ 06_svg/          # spatially variable genes

            âââ 07_communication/# CellChat results

            âââ 08_neighborhood/ # neighborhood enrichment



    Args:

        project_dir: Root path of the project to scaffold.

        overwrite: If False and the directory already exists, do nothing

            (returns a message). If True, recreate subdirectories.



    Returns:

        str: Log describing the created directory tree.

    """

    subdirs = [

        "raw_data/visium",

        "raw_data/scRNA",

        "raw_data/visium_hd",

        "raw_data/xenium",

        "raw_data/atera",

        "results/01_loading",

        "results/02_qc",

        "results/03_normalization",

        "results/04_clustering",

        "results/05_deconvolution",

        "results/06_svg",

        "results/07b_segmentation",

        "results/07_communication",

        "results/08_neighborhood",

        "results/10_velocity",

    ]

    if os.path.exists(project_dir) and not overwrite:

        return f"Project directory already exists: {project_dir} (use overwrite=True to recreate)"



    os.makedirs(project_dir, exist_ok=True)

    created = []

    for sub in subdirs:

        p = os.path.join(project_dir, sub)

        os.makedirs(p, exist_ok=True)

        created.append(p)



    readme = os.path.join(project_dir, "README.md")

    with open(readme, "w") as f:

        f.write("# Spatial transcriptomics project\n\n")

        f.write("## Directory layout\n\n")

        for sub in subdirs:

            f.write(f"- `{sub}/`\n")



    lines = [f"Initialized spatial project at {project_dir}"]

    lines.append("Created directories:")

    for p in created:

        lines.append(f"  {p}")

    lines.append(f"README written to {readme}")

    return "\n".join(lines)


# ============================================================================
# Visium HD / Xenium / Atera 平台加载与细胞重建
# ============================================================================

