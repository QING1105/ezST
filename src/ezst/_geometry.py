"""Data-only StarDist geometry and instance-label validation."""
from pathlib import Path
import numpy as np


def load_nuclei_polygons(path, *, allow_legacy_pickle=False):
    """Read numeric NPZ geometry; legacy pickle requires an explicit opt-in."""
    path = Path(path)
    if path.suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            polygons = {key: archive[key] for key in archive.files}
    elif allow_legacy_pickle and path.suffix.lower() in {".pckl", ".pkl", ".pickle"}:
        import pickle
        with path.open("rb") as handle:
            polygons = pickle.load(handle)
    else:
        raise ValueError(
            "Use numeric .npz geometry. For trusted legacy pickle files, "
            "pass allow_legacy_pickle=True explicitly."
        )
    if not isinstance(polygons, dict) or "coord" not in polygons:
        raise ValueError("StarDist geometry must contain a 'coord' array")
    coord = np.asarray(polygons["coord"])
    if coord.ndim != 3 or coord.shape[1] != 2 or coord.shape[2] < 3:
        raise ValueError("StarDist coord must have shape (n_nuclei, 2, n_vertices>=3)")
    if not np.issubdtype(coord.dtype, np.number) or not np.isfinite(coord).all():
        raise ValueError("StarDist coordinates must be finite numeric values")
    return polygons


def save_nuclei_polygons(path, polygons):
    path = Path(path)
    if path.suffix.lower() != ".npz":
        raise ValueError("Write StarDist geometry to a .npz file")
    arrays = {key: np.asarray(value) for key, value in polygons.items()}
    if any(array.dtype.hasobject for array in arrays.values()):
        raise ValueError("StarDist geometry must contain data-only arrays")
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **arrays)


def filter_nuclei_instances(labels, polygons, min_area):
    """Filter by original instance area, preserving touching nuclei as distinct."""
    labels = np.asarray(labels)
    if labels.ndim != 2 or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("Nucleus labels must be a two-dimensional integer array")
    if (labels < 0).any() or min_area < 0:
        raise ValueError("Nucleus labels and min_area must be nonnegative")
    n_polygons = len(polygons["coord"])
    if labels.max(initial=0) > n_polygons:
        raise ValueError("Nucleus label IDs exceed the number of polygons")
    sizes = np.bincount(labels.ravel(), minlength=n_polygons + 1)
    keep = np.flatnonzero((sizes[1:] > 0) & (sizes[1:] >= min_area))
    mapping = np.zeros(n_polygons + 1, dtype=np.uint32)
    mapping[keep + 1] = np.arange(1, len(keep) + 1, dtype=np.uint32)
    filtered = {}
    for key, value in polygons.items():
        array = np.asarray(value)
        filtered[key] = array[keep] if array.ndim and len(array) == n_polygons else array.copy()
    filtered["label_ids"] = np.arange(1, len(keep) + 1, dtype=np.uint32)
    return mapping[labels], filtered


def nuclei_xy_polygons(path, *, allow_legacy_pickle=False, object_ids=None):
    """Convert StarDist (row, col) vertices to plotting (x, y) coordinates."""
    polygons = load_nuclei_polygons(path, allow_legacy_pickle=allow_legacy_pickle)
    coord = polygons["coord"]
    if object_ids is not None:
        ids = np.asarray(object_ids, dtype=np.int64)
        if ((ids < 1) | (ids > len(coord))).any():
            raise ValueError("Cell object_id does not match the saved nucleus polygons")
        coord = coord[ids - 1]
    return [vertices[::-1].T for vertices in coord]


def read_proseg_counts(counts_path, cells_path, genes_path, *, microns_per_pixel):
    """Read Proseg's cell-by-gene matrix and its exported metadata in order."""
    import anndata as ad
    import pandas as pd
    from scipy.io import mmread
    if not np.isfinite(microns_per_pixel) or microns_per_pixel <= 0:
        raise ValueError("microns_per_pixel must be finite and positive")
    counts = mmread(counts_path).tocsr()
    cells = pd.read_csv(cells_path)
    genes = pd.read_csv(genes_path)
    required = {"cell", "centroid_x", "centroid_y"}
    if not required.issubset(cells.columns) or "gene" not in genes.columns:
        raise ValueError("Proseg metadata requires cell/centroid_x/centroid_y and gene columns")
    if counts.shape != (len(cells), len(genes)):
        raise ValueError(
            f"Proseg matrix {counts.shape} does not match {len(cells)} cells and {len(genes)} genes"
        )
    cells = cells.sort_values("cell").reset_index(drop=True)
    if not np.array_equal(cells["cell"].to_numpy(), np.arange(len(cells))):
        raise ValueError("Proseg cell IDs must identify matrix rows 0..n_cells-1 exactly once")
    if genes["gene"].isna().any() or genes["gene"].astype(str).str.strip().eq("").any():
        raise ValueError("Proseg gene metadata contains missing gene names")
    if genes["gene"].duplicated().any():
        raise ValueError("Proseg gene metadata contains duplicate gene names")
    coords = cells[["centroid_x", "centroid_y"]].to_numpy(dtype=float)
    if not np.isfinite(coords).all():
        raise ValueError("Proseg cell centroids must be finite")
    if not np.isfinite(counts.data).all() or (counts.data < 0).any():
        raise ValueError("Proseg counts must be finite and nonnegative")
    cells.index = pd.Index([f"cell_{cell}" for cell in cells["cell"]], name="barcode")
    genes.index = pd.Index(genes["gene"].astype(str), name="gene_name")
    result = ad.AnnData(X=counts, obs=cells, var=genes)
    result.obsm["spatial"] = coords / microns_per_pixel
    result.uns["cellpose_scale"] = float(microns_per_pixel)
    result.uns["proseg"] = {"counts_path": str(counts_path), "genes_path": str(genes_path),
                            "matrix_order": "cells_by_genes", "coordinate_units": "fullres_pixels"}
    return result
