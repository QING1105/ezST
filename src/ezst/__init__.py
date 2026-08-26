"""
ezST — Independent spatial transcriptomics analysis toolkit.

Self-contained implementation of the classic 10x Visium (spot-level) workflow:
load -> QC -> normalize -> cluster -> domains/SVG -> deconvolution (DestVI/SPOTlight)
-> neighborhood enrichment / CellChat.

No external runtime dependency (does NOT require Biomni).

Platform extensions (Visium HD / Xenium / Atera) are under active development
in `ezst.visium_new_platforms`.
"""

__version__ = "0.2.0"

from .visium import (  # noqa: F401
    init_spatial_project,
    load_visium_data,
    filter_visium_spots,
    normalize_visium,
    cluster_spatial_data,
    identify_spatial_domains,
    find_spatially_variable_genes,
    deconvolve_spatial_destvi,
    deconvolve_spatial_spotlight,
    spatial_neighborhood_enrichment,
    infer_spatial_cell_communication,
)
from .visium_new_platforms import (  # noqa: F401
    load_visium_hd,
    segment_cells,
    bin_to_cell,
    run_proseg,
    load_xenium,
    load_atera,
    plot_cells_crop,
)
