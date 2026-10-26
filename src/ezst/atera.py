"""Read official Atera AnnData Zarr matrices without loading all counts into RAM."""
from contextlib import contextmanager
from pathlib import Path
import json


@contextmanager
def open_atera_zarr(path):
    import zarr
    path = Path(path)
    if path.suffix == ".zip":
        with zarr.storage.ZipStore(str(path), mode="r") as store:
            yield zarr.open_group(store=store, mode="r")
    else:
        yield zarr.open_group(str(path), mode="r")


def _resolve_input(source):
    source = Path(source).resolve()
    directory = source if source.is_dir() and not source.name.endswith('.zarr') else source.parent
    manifests = list(directory.glob('*experiment.spatial'))
    if source.name.endswith('experiment.spatial'):
        manifest_path = source
    elif len(manifests) == 1:
        manifest_path = manifests[0]
    elif len(manifests) > 1:
        raise ValueError('Multiple experiments found; pass one experiment.spatial path')
    else:
        manifest_path = None
    manifest = json.loads(manifest_path.read_text(encoding='utf-8')) if manifest_path else {}
    if (source.is_file() and source.suffix == '.zip') or source.name.endswith('.zarr'):
        matrix = source
    else:
        relative = manifest.get('explorer_files', {}).get('cell_feature_zarr_filepath', 'cell_feature_matrix.zarr.zip')
        matrix = directory / relative
        if not matrix.exists():
            candidates = [p for p in directory.glob('*cell_feature_matrix.zarr.zip')
                          if not p.name.endswith('csc_cell_feature_matrix.zarr.zip')]
            if len(candidates) != 1:
                raise FileNotFoundError('Supply cell_feature_matrix.zarr.zip and its experiment.spatial')
            matrix = candidates[0]
    if matrix.name.endswith('csc_cell_feature_matrix.zarr.zip'):
        raise ValueError('Use cell_feature_matrix.zarr.zip, which includes cell metadata and analysis results')
    return matrix, manifest


def _read_elem(group):
    try:
        from anndata.io import read_elem
    except ImportError:
        from anndata._io.specs import read_elem
    return read_elem(group)


def load_atera_zarr(source, output_path, sample_id=None, roi=None,
                    use_official_filter=True, cell_type_key='azimuth_coarse',
                    chunk_rows=2048):
    """Stream genes into h5ad; ROI is (xmin, ymin, xmax, ymax), in micrometres.

    All source rows are checked against official gene/control totals. The output
    contains the selected cells and Gene Expression columns; control summaries,
    original annotations and feature IDs are retained. X and layers/counts are
    HDF5 hard links on disk, saving a second copy of the original counts.
    """
    import anndata as ad
    import h5py
    import numpy as np
    import pandas as pd
    from scipy import sparse

    if not isinstance(chunk_rows, int) or chunk_rows < 1:
        raise ValueError('chunk_rows must be a positive integer')
    matrix, manifest = _resolve_input(source)
    output = Path(output_path).resolve()
    if output == matrix:
        raise ValueError('Output must differ from input')
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix('.partial.h5ad')
    with open_atera_zarr(matrix) as group:
        if group['X'].attrs.get('encoding-type') != 'csr_matrix':
            raise ValueError('Atera cell_feature_matrix must contain a CSR matrix')
        shape = tuple(group['X'].attrs['shape'])
        obs, var = _read_elem(group['obs']), _read_elem(group['var'])
        if shape != (len(obs), len(var)) or not obs.index.is_unique:
            raise ValueError('Matrix dimensions or cell barcodes are inconsistent')
        coords = obs[['centroid_column', 'centroid_row']].to_numpy(dtype=float)
        if not np.isfinite(coords).all():
            raise ValueError('Cell centroids must be finite micrometre coordinates')
        gene_mask = var['feature_type'].astype(str).eq('Gene Expression').to_numpy()
        if not gene_mask.any():
            raise ValueError('No Gene Expression features found')
        keep = np.ones(len(obs), dtype=bool)
        if use_official_filter:
            if 'filtered' not in obs or not pd.api.types.is_bool_dtype(obs['filtered']):
                raise ValueError('Official filtered field must be boolean (True means retained)')
            keep &= obs['filtered'].to_numpy()
        if roi is not None:
            roi = np.asarray(roi, dtype=float)
            if roi.shape != (4,) or not np.isfinite(roi).all() or np.any(roi[2:] <= roi[:2]):
                raise ValueError('roi must be xmin, ymin, xmax, ymax in micrometres')
            keep &= ((coords >= roi[:2]) & (coords < roi[2:])).all(axis=1)
        if not keep.any():
            raise ValueError('The selected ROI/filter contains no cells')
        selected = obs.loc[keep].copy()
        selected['official_keep'] = obs.loc[keep, 'filtered'].astype(bool)
        selected['sample_id'] = sample_id or manifest.get('run_name', matrix.parent.name)
        selected['roi_id'] = manifest.get('roi_uuid', manifest.get('region_name', 'roi_0'))
        if cell_type_key is not None:
            if cell_type_key not in selected:
                raise ValueError(f'Cell type field {cell_type_key!r} is absent; pass cell_type_key=None to load without labels')
            selected['cell_type'] = pd.Categorical(selected[cell_type_key].astype(str))
        genes = var.loc[gene_mask].copy()
        genes.index = pd.Index(genes['feature_name'].astype(str), name='gene_symbol')
        meta = {
            'platform': 'atera', 'coordinate_unit': 'micrometers',
            'source': str(matrix), 'source_n_cells': len(obs), 'source_n_features': len(var),
            'use_official_filter': bool(use_official_filter),
            'cell_type_source': cell_type_key or '',
            'experiment_json': json.dumps(manifest, ensure_ascii=False),
            'official_zarr_metadata_json': json.dumps(dict(group.attrs), ensure_ascii=False),
        }
        if roi is not None:
            meta['roi_bounds_um'] = roi
        embeddings = {'spatial': coords[keep]}
        if 'obsm' in group:
            for key in group['obsm'].keys():
                embedding = np.asarray(group['obsm'][key])
                if len(embedding) == len(obs):
                    embeddings[key + '_10x'] = embedding[keep]
        result = ad.AnnData(obs=selected, var=genes, obsm=embeddings, uns={'spatial_meta': meta})
        result.var_names_make_unique()
        result.write_h5ad(temporary)
        del result
        pointers = np.asarray(group['X/indptr'], dtype=np.int64)
        gene_counts = np.zeros(len(obs), dtype=np.int64)
        n_genes = np.zeros(len(obs), dtype=np.int32)
        total_counts = np.zeros(len(obs), dtype=np.int64)
        written_rows = written_nnz = 0
        with h5py.File(temporary, 'a') as handle:
            x = handle.create_group('X')
            x.attrs.update({'encoding-type': 'csr_matrix', 'encoding-version': '0.1.0',
                            'shape': np.asarray([int(keep.sum()), int(gene_mask.sum())])})
            values = x.create_dataset('data', shape=(0,), maxshape=(None,), dtype='int32', compression='lzf')
            indices = x.create_dataset('indices', shape=(0,), maxshape=(None,), dtype='int32', compression='lzf')
            indptr = x.create_dataset('indptr', shape=(int(keep.sum())+1,), dtype='int64', compression='lzf')
            indptr[0] = 0
            for start in range(0, len(obs), chunk_rows):
                end = min(start + chunk_rows, len(obs))
                lo, hi = int(pointers[start]), int(pointers[end])
                data = np.asarray(group['X/data'][lo:hi])
                cols = np.asarray(group['X/indices'][lo:hi])
                if data.dtype.kind not in 'iu' or np.any(data < 0) or np.any(data > np.iinfo(np.int32).max):
                    raise ValueError('Atera counts must be nonnegative int32 values')
                block = sparse.csr_matrix((data, cols, pointers[start:end+1]-lo), shape=(end-start, len(var)))
                total_counts[start:end] = np.asarray(block.sum(axis=1)).ravel()
                block = block[:, gene_mask]
                block.eliminate_zeros()
                gene_counts[start:end] = np.asarray(block.sum(axis=1)).ravel()
                n_genes[start:end] = block.getnnz(axis=1)
                block = block[keep[start:end]]
                if not block.shape[0]:
                    continue
                next_nnz = written_nnz + block.nnz
                values.resize((next_nnz,)); indices.resize((next_nnz,))
                values[written_nnz:next_nnz] = block.data
                indices[written_nnz:next_nnz] = block.indices
                indptr[written_rows+1:written_rows+block.shape[0]+1] = block.indptr[1:] + written_nnz
                written_rows += block.shape[0]
                written_nnz = next_nnz
            for key, actual in [('transcript_counts', gene_counts), ('total_counts', total_counts)]:
                if key in obs and not np.array_equal(actual, obs[key].to_numpy(dtype=np.int64)):
                    raise ValueError(f'Matrix counts disagree with official obs[{key!r}]')
            handle['layers']['counts'] = x
            selected['n_genes_by_counts'] = n_genes[keep]
            selected['total_gene_counts'] = gene_counts[keep]
            del handle['obs']
            try:
                from anndata.io import write_elem
            except ImportError:
                from anndata._io.specs import write_elem
            write_elem(handle, 'obs', selected)
        temporary.replace(output)
        summary = {
            'source_cells': len(obs), 'source_features': len(var), 'gene_features': int(gene_mask.sum()),
            'control_features': int((~gene_mask).sum()),
            'official_retained_cells': int(obs['filtered'].sum()), 'output_cells': int(keep.sum()),
            'source_nonzero_entries': int(pointers[-1]), 'output_gene_nonzero_entries': written_nnz,
            'gene_counts_sum': int(gene_counts.sum()), 'all_feature_counts_sum': int(total_counts.sum()),
            'negative_control_probe_counts': int(obs['control_probe_counts'].sum()),
            'negative_control_codeword_counts': int(obs['control_codeword_counts'].sum()),
            'median_genes_positive_cells': float(np.median(n_genes[gene_counts > 0])),
            'median_transcripts_positive_cells': float(np.median(gene_counts[gene_counts > 0])),
            'source_count_checks_passed': True, 'coordinate_unit': 'micrometers',
            'roi_bounds_um': roi.tolist() if roi is not None else None,
            'analysis_version': manifest.get('analysis_sw_version', ''),
            'output': str(output),
        }
        output.with_suffix('.qc.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
        qc = obs[['centroid_column', 'centroid_row', 'transcript_counts', 'total_counts',
                  'control_probe_counts', 'control_codeword_counts', 'cell_area', 'nucleus_area', 'filtered']].copy()
        qc['n_genes_by_counts'] = n_genes
        qc.to_csv(output.with_suffix('.source_qc.csv.gz'), compression='gzip')
    return json.dumps(summary, ensure_ascii=False, indent=2)


def _plot_metadata(adata_path, embeddings=('spatial',)):
    """Read plot metadata only: AnnData backed mode still eagerly reads layers."""
    import anndata as ad
    import h5py
    with h5py.File(adata_path, 'r') as handle:
        return ad.AnnData(obs=_read_elem(handle['obs']),
                         obsm={key: _read_elem(handle['obsm'][key]) for key in embeddings},
                         uns={'spatial_meta': _read_elem(handle['uns/spatial_meta'])})


def plot_atera_qc(adata_path, output_dir):
    """Plot full-source count QC and the actual selected ROI using physical units."""
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    a = _plot_metadata(adata_path)
    qc = pd.read_csv(Path(adata_path).with_suffix('.source_qc.csv.gz'), index_col=0)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, column, title in zip(axes[0], ['transcript_counts', 'n_genes_by_counts', 'cell_area'],
                                 ['Gene transcripts / cell', 'Detected genes / cell', 'Cell area (um2)']):
        ax.hist(qc[column], bins=80, color='#3287a8'); ax.set(title=title, ylabel='Cells')
    view = qc.iloc[::max(1, len(qc)//80000)]
    p = axes[1, 0].scatter(view.centroid_column, view.centroid_row, s=.2,
                           c=np.log1p(view.transcript_counts), cmap='viridis', rasterized=True)
    fig.colorbar(p, ax=axes[1, 0], label='log1p gene transcripts')
    axes[1, 0].set_title(f'Full tissue: {len(qc):,} cells')
    if 'roi_bounds_um' in a.uns['spatial_meta']:
        from matplotlib.patches import Rectangle
        x0, y0, x1, y1 = a.uns['spatial_meta']['roi_bounds_um']
        axes[1, 0].add_patch(Rectangle((x0, y0), x1-x0, y1-y0, fill=False,
                                      edgecolor='#d55e00', linewidth=1.5))
    low = qc.loc[~qc.filtered.astype(bool)]
    axes[1, 1].scatter(view.centroid_column, view.centroid_row, s=.1, c='#cccccc', rasterized=True)
    axes[1, 1].scatter(low.centroid_column, low.centroid_row, s=4, c='#d55e00', label=f'Official excluded: {len(low)}')
    axes[1, 1].legend(fontsize=8); axes[1, 1].set_title('Official cell filter')
    coords = a.obsm['spatial']
    axes[1, 2].scatter(coords[:, 0], coords[:, 1], s=2, c=a.obs['total_gene_counts'], cmap='viridis', rasterized=True)
    axes[1, 2].set_title(f'Selected analysis region: {a.n_obs:,} cells')
    for ax in axes[1]:
        ax.set(xlabel='x (um)', ylabel='y (um)', aspect='equal'); ax.invert_yaxis()
    fig.tight_layout(); fig.savefig(output/'atera_qc.png', dpi=160); plt.close(fig)
    return str(output/'atera_qc.png')


def plot_atera_segmentation(adata_path, cells_path, output_path):
    """Match official segmentation by cell ID and plot polygons in micrometres."""
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    a = _plot_metadata(adata_path)
    with open_atera_zarr(cells_path) as cells:
        ids = np.asarray(cells['cell_id'], dtype=np.uint64)
        combined = ids[:, 0] | (ids[:, 1] << np.uint64(32))
        target = a.obs['cell_id'].to_numpy(dtype=np.uint64)
        order = np.argsort(combined)
        found = np.searchsorted(combined[order], target)
        if np.any(found >= len(combined)) or not np.array_equal(combined[order[found]], target):
            raise ValueError('Matrix and segmentation cell IDs do not match')
        cell_indices = order[found]
        columns = list(cells['cell_summary'].attrs['column_names'])
        summaries = np.asarray(cells['cell_summary'])[cell_indices]
        centers = summaries[:, [columns.index('cell_centroid_x'), columns.index('cell_centroid_y')]]
        error = float(np.max(np.abs(centers - a.obsm['spatial'])))
        if error > .01:
            raise ValueError(f'Matrix/segmentation coordinates differ by {error} micrometres')
        polygon_set = cells['polygon_sets/1']
        polygon_ids = np.asarray(polygon_set['cell_index'])
        rows = np.flatnonzero(np.isin(polygon_ids, cell_indices))
        n_vertices = np.asarray(polygon_set['num_vertices'])[rows]
        vertices = np.asarray(polygon_set['vertices'].oindex[rows, :])
        polygons = [v[:2*int(n)].reshape(-1, 2) for v, n in zip(vertices, n_vertices) if n >= 3]
        transform = np.asarray(cells['masks/homogeneous_transform'])
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.add_collection(PolyCollection(polygons, facecolors='#dbeaf0', edgecolors='#285c75', linewidths=.35))
    ax.scatter(centers[:, 0], centers[:, 1], c='#d55e00', s=.5, label='Matrix cell centroids')
    ax.autoscale(); ax.set_aspect('equal'); ax.invert_yaxis()
    ax.set(xlabel='x (um)', ylabel='y (um)', title=f'Atera official segmentation: {a.n_obs:,} matched cells')
    ax.legend(loc='upper right', fontsize=8)
    fig.tight_layout(); fig.savefig(output, dpi=180); plt.close(fig)
    report = {'matched_cells': a.n_obs, 'plotted_polygons': len(polygons),
              'max_centroid_error_um': error, 'microns_to_pixels': transform.tolist()}
    output.with_suffix('.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def _annotation_geometry(a, cells_path):
    """Return polygons and their observation rows; align by ID, never by row order."""
    import numpy as np
    with open_atera_zarr(cells_path) as cells:
        ids = np.asarray(cells['cell_id'], dtype=np.uint64)
        combined = ids[:, 0] | (ids[:, 1] << np.uint64(32))
        target = a.obs['cell_id'].to_numpy(dtype=np.uint64)
        if len(np.unique(combined)) != len(combined) or len(np.unique(target)) != len(target):
            raise ValueError('Cell IDs must be unique in matrix and segmentation')
        order = np.argsort(combined)
        found = np.searchsorted(combined[order], target)
        if np.any(found >= len(combined)) or not np.array_equal(combined[order[found]], target):
            raise ValueError('Matrix and segmentation cell IDs do not match')
        matched = order[found]
        columns = list(cells['cell_summary'].attrs['column_names'])
        summary = np.asarray(cells['cell_summary'])[matched]
        centers = summary[:, [columns.index('cell_centroid_x'), columns.index('cell_centroid_y')]]
        error = float(np.max(np.abs(centers - a.obsm['spatial'])))
        if error > .01:
            raise ValueError(f'Matrix/segmentation coordinates differ by {error} micrometres')
        obs_for_cell = np.full(len(combined), -1, dtype=np.int64)
        obs_for_cell[matched] = np.arange(a.n_obs)
        group = cells['polygon_sets/1']
        polygon_cells = np.asarray(group['cell_index'], dtype=np.int64)
        sizes = np.asarray(group['num_vertices'])
        rows = np.flatnonzero((obs_for_cell[polygon_cells] >= 0) & (sizes >= 3))
        vertices = np.asarray(group['vertices'].oindex[rows, :])
        polygons = [v[:2*int(n)].reshape(-1, 2) for v, n in zip(vertices, sizes[rows])]
        obs_rows = obs_for_cell[polygon_cells[rows]]
    return polygons, obs_rows, rows, error


def _annotation_palette(labels):
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_hex
    palette = {
        'Renal epithelial cell - proximal tubules': '#E9A31B',
        'Renal epithelial cell - Loop of Henle': '#397AC0',
        'Renal epithelial cell - distal tubules': '#8752AB',
        'Podocyte': '#C33163', 'Type A intercalated cell': '#AD550B',
        'Type B intercalated cell': '#E383B2', 'Endothelial cell': '#DE4039',
        'Pericyte': '#33A5AF', 'Vascular smooth muscle cell': '#12696B',
        'Fibroblast': '#8C7353', 'Macrophage': '#8DA927', 'T cell': '#087C4F',
        'B cell': '#58BD8C', 'Dendritic cell': '#4B5731', 'Granulocyte': '#C6CA3A',
        'Red blood cell': '#6D283D', 'Coarse Not Confidently Annotated': '#BFBFBF',
        'Unassigned': '#7E7E7E', 'Unannotated': '#7E7E7E',
    }
    extras = sorted(set(labels) - palette.keys())
    for i, label in enumerate(extras):
        palette[label] = to_hex(plt.get_cmap('turbo')((i+.5)/max(1, len(extras))))
    return {label: palette[label] for label in sorted(set(labels))}


def _short_cell_label(label):
    return (label.replace('Renal epithelial cell - ', '').replace('Epithelial cell of ', 'Epithelial: ')
            .replace('Coarse Not Confidently Annotated', 'Low-confidence / coarse')
            .replace('Vascular smooth muscle cell', 'Vascular smooth muscle')
            .replace('Hematopoietic precursor cell', 'Hematopoietic precursor')
            .replace('Type A intercalated cell', 'Intercalated A')
            .replace('Type B intercalated cell', 'Intercalated B'))


def plot_atera_annotations(adata_path, output_path, cells_path=None, source_matrix_path=None,
                           label_key='azimuth_coarse', title='Atera sample'):
    """Draw spatial annotation on cell polygons, plus a complete tissue map.

    Both plots share one label-to-color mapping. Sidecars retain the original
    labels, cell IDs, polygon rows and colors so the drawing is auditable.
    """
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Rectangle
    a = _plot_metadata(adata_path)
    source = Path(source_matrix_path or a.uns['spatial_meta']['source'])
    cells_path = Path(cells_path) if cells_path else source.parent/'cells.zarr.zip'
    with open_atera_zarr(source) as g:
        full = _read_elem(g['obs'])
    full = full.loc[full['filtered']].copy()
    labels = a.obs[label_key].astype('string').fillna('Unassigned')
    full_labels = full[label_key].astype('string').fillna('Unassigned')
    full_by_id = pd.Series(full_labels.to_numpy(), index=full['cell_id'].to_numpy())
    if not np.array_equal(full_by_id.loc[a.obs['cell_id']].to_numpy(), labels.to_numpy()):
        raise ValueError('ROI annotations differ from the selected source annotation field')
    polygons, obs_rows, polygon_rows, error = _annotation_geometry(a, cells_path)
    palette = _annotation_palette(set(labels) | set(full_labels))
    colors = np.asarray([palette[label] for label in labels])
    polygon_colors = colors[obs_rows]
    counts = labels.value_counts()
    unresolved = [s for s in counts.index if s in {'Unassigned', 'Unannotated'} or 'Not Confidently Annotated' in s]
    categories = [s for s in counts.index if s not in unresolved] + unresolved
    xy = a.obsm['spatial']
    bounds = np.asarray(a.uns['spatial_meta'].get('roi_bounds_um', [*xy.min(0), *xy.max(0)]))
    x0, y0, x1, y1 = bounds
    center = (bounds[:2] + bounds[2:])/2
    half = (bounds[2:] - bounds[:2])/6
    detail = np.r_[center-half, center+half]
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    def spatial_axes(ax, box, ticks=True):
        ax.set(xlim=(box[0], box[2]), ylim=(box[3], box[1]), aspect='equal')
        ax.set_xlabel('x (um)', fontsize=10); ax.set_ylabel('y (um)', fontsize=10)
        ax.tick_params(labelsize=9)
        if not ticks:
            ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(''); ax.set_ylabel('')

    def draw_cells(ax, box):
        ax.add_collection(PolyCollection(polygons, facecolors=polygon_colors,
                                        edgecolors='#ffffff', linewidths=.16, rasterized=True))
        spatial_axes(ax, box)

    def legend(ax, order, frequencies, columns=1, fontsize=10):
        ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        per_col = int(np.ceil(len(order)/columns))
        step = min(.034, .86/max(1, per_col))
        for col in range(columns):
            left, right = col/columns, (col+1)/columns
            ax.text(left+.035, .975, 'Cell type', fontsize=fontsize+1, weight='bold')
            ax.text(right-.025, .975, 'Cells', fontsize=fontsize+1, weight='bold', ha='right')
        for i, label in enumerate(order):
            col, row = divmod(i, per_col)
            left, right = col/columns, (col+1)/columns
            y = .94-row*step
            ax.add_patch(Rectangle((left+.008, y-.009), .018, .018,
                                   facecolor=palette[label], edgecolor='none'))
            ax.text(left+.036, y, _short_cell_label(label), fontsize=fontsize, va='center')
            ax.text(right-.025, y, f'{int(frequencies[label]):,}', fontsize=fontsize,
                    va='center', ha='right')

    fig = plt.figure(figsize=(21, 11.5), facecolor='white')
    grid = fig.add_gridspec(2, 3, width_ratios=[.65, 1.8, 1.04], height_ratios=[.8, 1.2],
                           left=.035, right=.99, top=.90, bottom=.075, wspace=.18, hspace=.2)
    locator = fig.add_subplot(grid[0, 0])
    location = full.iloc[::max(1, len(full)//100000)]
    locator.scatter(location.centroid_column, location.centroid_row, s=.2, color='#bdc5cb', rasterized=True)
    locator.add_patch(Rectangle((x0, y0), x1-x0, y1-y0, fill=False, edgecolor='#111111', linewidth=1.8))
    locator.set_title('A  Region within the sample', fontsize=12, loc='left', pad=12)
    locator.set_aspect('equal'); locator.invert_yaxis(); locator.set_xticks([]); locator.set_yticks([])
    detail_ax = fig.add_subplot(grid[1, 0]); draw_cells(detail_ax, detail)
    detail_ax.set_title('C  Cell-boundary detail', fontsize=12, loc='left', pad=12)
    main = fig.add_subplot(grid[:, 1]); draw_cells(main, bounds)
    main.add_patch(Rectangle(detail[:2], detail[2]-detail[0], detail[3]-detail[1],
                             fill=False, edgecolor='#333333', linewidth=1, linestyle='--'))
    main.set_title(f'B  Cell types in tissue space | {a.n_obs:,} cells', fontsize=16, loc='left', pad=13)
    legend_ax = fig.add_subplot(grid[:, 2]); legend(legend_ax, categories, counts)
    legend_ax.text(.035, .035, f'{len(categories)} original labels | {a.n_obs:,} cells\nGray: unresolved official annotations',
                   fontsize=10, color='#555555', va='bottom')
    fig.suptitle(f'{title} | Spatial cell-type annotation', x=.035, ha='left', fontsize=23, weight='bold', y=.972)
    fig.text(.035, .932, f'Official {label_key} labels on 10x cell boundaries; colors identify cell types.', fontsize=12, color='#444444')
    fig.text(.035, .025, f'ROI: x = {x0:g}-{x1:g} um, y = {y0:g}-{y1:g} um  |  Cell IDs and centroids matched against cells.zarr.zip',
             fontsize=11, color='#444444')
    fig.savefig(output, dpi=190); plt.close(fig)

    full_output = output.with_name(output.stem+'_full_tissue.png')
    full_counts = full_labels.value_counts()
    full_order = [s for s in full_counts.index if s not in unresolved] + [s for s in unresolved if s in full_counts]
    fig = plt.figure(figsize=(23, 11), facecolor='white')
    ax = fig.add_axes([.035, .14, .54, .74])
    ax.scatter(full.centroid_column, full.centroid_row,
               c=[palette[label] for label in full_labels], s=.45, linewidths=0, rasterized=True)
    ax.add_patch(Rectangle((x0, y0), x1-x0, y1-y0, fill=False, edgecolor='#111111', linewidth=1.5))
    ax.annotate('Detailed ROI', (x1, y0), xytext=(30, 30), textcoords='offset points', fontsize=11,
                arrowprops={'arrowstyle':'->', 'color':'#222222'})
    ax.set_aspect('equal'); ax.invert_yaxis(); ax.set(xlabel='x (um)', ylabel='y (um)')
    legend_ax = fig.add_axes([.60, .10, .39, .80])
    legend(legend_ax, full_order, full_counts, columns=2, fontsize=8.5)
    fig.suptitle(f'{title} | Full-tissue spatial annotation', x=.035, ha='left', fontsize=23, weight='bold', y=.965)
    fig.text(.035, .915, f'{len(full):,} officially retained cells | {len(full_counts)} original {label_key} labels | same colors as the ROI', fontsize=12)
    fig.text(.035, .045, 'Each point is one cell at its measured tissue coordinate. The black box locates the detailed cell-boundary view.', fontsize=11, color='#444444')
    fig.savefig(full_output, dpi=190); plt.close(fig)

    mapping = pd.DataFrame({'barcode': a.obs_names, 'cell_id': a.obs['cell_id'].to_numpy(),
                            'x_um': xy[:, 0], 'y_um': xy[:, 1], 'official_label': labels.to_numpy(),
                            'color': colors, 'polygon_row': -1})
    mapping.loc[obs_rows, 'polygon_row'] = polygon_rows
    mapping.to_csv(output.with_suffix('.cells.csv'), index=False)
    pd.DataFrame({'official_label': full_counts.index, 'display_label': [_short_cell_label(s) for s in full_counts.index],
                  'cells': full_counts.to_numpy(), 'color': [palette[s] for s in full_counts.index]}).to_csv(
                      output.with_suffix('.legend.csv'), index=False)
    report = {'annotation_key': label_key, 'source_cells': len(full), 'roi_cells': a.n_obs,
              'roi_label_counts': {str(k): int(v) for k, v in counts.items()},
              'full_label_count': len(full_counts), 'plotted_polygons': len(polygons),
              'polygon_cell_ids_unique': bool(mapping.loc[mapping.polygon_row >= 0, 'cell_id'].is_unique),
              'max_centroid_error_um': error, 'roi_bounds_um': bounds.tolist(),
              'roi_figure': str(output), 'full_tissue_figure': str(full_output), 'palette': palette}
    output.with_suffix('.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return str(output)
