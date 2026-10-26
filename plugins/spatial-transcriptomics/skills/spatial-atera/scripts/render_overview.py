"""Render the kidney example's completed stages as one scientific figure."""
from pathlib import Path
import json
import sys


def render_overview(data, run, output):
    import numpy as np
    import pandas as pd
    import anndata as ad
    import h5py
    import networkx as nx
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Rectangle
    from matplotlib.colors import Normalize
    try:
        from ezst.atera import _plot_metadata, _annotation_geometry, _short_cell_label, open_atera_zarr, _read_elem
    except ImportError:
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'spatial-visium-hd'))
        from atera import _plot_metadata, _annotation_geometry, _short_cell_label, open_atera_zarr, _read_elem

    data, run, output = Path(data), Path(run), Path(output)
    loaded = _plot_metadata(run/'01_loading/kidney_cells.h5ad')
    clustered = _plot_metadata(run/'02_expression/kidney_clustered.h5ad', ('spatial', 'X_umap'))
    domains = _plot_metadata(run/'03_domains/kidney_domains.h5ad')
    for result in (clustered, domains):
        if not result.obs_names.equals(loaded.obs_names):
            raise ValueError('Overview stages must contain the same cells in the same order')
        np.testing.assert_array_equal(result.obsm['spatial'], loaded.obsm['spatial'])
    qc = pd.read_csv(run/'01_loading/kidney_cells.source_qc.csv.gz', index_col=0)
    svg = pd.read_csv(run/'03_domains/svg.csv', index_col=0)
    enrichment = pd.read_csv(run/'04_neighborhood/enrichment.csv', index_col=0)
    communication = pd.read_csv(run/'05_communication/cellchat_communication.csv')
    selection = json.loads((run/'04_neighborhood/cell_selection.json').read_text())
    annotation = json.loads((run/'02_expression/official_annotations.json').read_text())
    palette = annotation['palette']
    with open_atera_zarr(data/'cell_feature_matrix.zarr.zip') as group:
        full = _read_elem(group['obs'])
    full = full.loc[full['filtered']]
    full_labels = full['azimuth_coarse'].astype('string').fillna('Unassigned')
    labels = loaded.obs['azimuth_coarse'].astype('string').fillna('Unassigned')
    counts = labels.value_counts()
    unresolved = [s for s in counts.index if 'Not Confidently Annotated' in s or s in {'Unassigned', 'Unannotated'}]
    order = [s for s in counts.index if s not in unresolved] + unresolved
    polygons, obs_rows, _, error = _annotation_geometry(loaded, data/'cells.zarr.zip')
    bounds = annotation['roi_bounds_um']
    x0, y0, x1, y1 = bounds
    abbreviations = {
        'Renal epithelial cell - proximal tubules': 'PT',
        'Renal epithelial cell - Loop of Henle': 'LH',
        'Renal epithelial cell - distal tubules': 'DT',
        'Endothelial cell': 'Endo', 'Macrophage': 'Mac', 'Pericyte': 'Peri',
        'T cell': 'T', 'Podocyte': 'Podo', 'Fibroblast': 'Fib',
        'Dendritic cell': 'DC', 'B cell': 'B', 'Vascular smooth muscle cell': 'VSM',
        'Type B intercalated cell': 'IC-B', 'Red blood cell': 'RBC',
    }
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':11,
                         'axes.titlesize':16, 'axes.titleweight':'bold',
                         'axes.labelsize':11, 'xtick.labelsize':9, 'ytick.labelsize':9})
    fig = plt.figure(figsize=(24, 21), facecolor='white')
    grid = fig.add_gridspec(3, 3, left=.045, right=.975, bottom=.085, top=.915,
                           height_ratios=[1.15, 1, 1], wspace=.28, hspace=.36)
    fig.suptitle('Atera pipeline | FFPE human kidney', x=.045, y=.975,
                 ha='left', fontsize=28, weight='bold')
    fig.text(.045, .949, 'Official 10x Atera 1.0.4 data  |  Cell types shown at their measured tissue positions', fontsize=15)
    fig.text(.045, .929, f'Full tissue: {len(full):,} cells  |  ROI: {loaded.n_obs:,} cells, x={x0:g}-{x1:g}, y={y0:g}-{y1:g} um', fontsize=12, color='#555555')

    def spatial(ax):
        ax.set(xlim=(x0, x1), ylim=(y1, y0), aspect='equal', xlabel='x (um)', ylabel='y (um)')

    def polygon_panel(ax, colors):
        ax.add_collection(PolyCollection(polygons, facecolors=np.asarray(colors)[obs_rows],
                                        edgecolors='white', linewidths=.12, rasterized=True))
        spatial(ax)

    ax = fig.add_subplot(grid[0, 0])
    ax.scatter(full.centroid_column, full.centroid_row, c=[palette[s] for s in full_labels],
               s=.25, linewidths=0, rasterized=True)
    ax.add_patch(Rectangle((x0, y0), x1-x0, y1-y0, fill=False, edgecolor='black', linewidth=1.8))
    ax.annotate('ROI in B', (x1, y0), xytext=(24, 24), textcoords='offset points',
                arrowprops={'arrowstyle':'->'}, fontsize=11)
    ax.set_aspect('equal'); ax.invert_yaxis()
    ax.set(xlabel='x (um)', ylabel='y (um)')
    ax.set_title('A  Full-tissue spatial cell types', loc='left', pad=13)
    ax.text(0, -.16, f'{full_labels.nunique()} original labels; same colors as B.\nFull legend accompanies the example.', transform=ax.transAxes,
            fontsize=11, color='#555555', va='top')

    ax = fig.add_subplot(grid[0, 1])
    polygon_panel(ax, [palette[s] for s in labels])
    ax.set_title('B  Cell types on ROI boundaries', loc='left', pad=13)

    ax = fig.add_subplot(grid[0, 2]); ax.set_axis_off()
    ax.set_title('C  ROI cell types and counts', loc='left', pad=13)
    for i, label in enumerate(order):
        y = .96-i*.037
        ax.add_patch(Rectangle((.01, y-.009), .022, .018, facecolor=palette[label], edgecolor='none'))
        text = _short_cell_label(label)
        if label in abbreviations:
            text += ' ('+abbreviations[label]+')'
        ax.text(.05, y, text, va='center', fontsize=11)
        ax.text(.97, y, f'{counts[label]:,}', ha='right', va='center', fontsize=11)
    ax.text(.01, -.025, 'Gray: low-confidence / unassigned official labels.\nAbbreviations are reused in H and I.', fontsize=10, color='#555555')

    ax = fig.add_subplot(grid[1, 0])
    cluster_labels = clustered.obs.leiden.astype('category')
    cluster_colors = plt.get_cmap('tab20')(np.arange(len(cluster_labels.cat.categories)))
    xy = clustered.obsm['X_umap']
    ax.scatter(xy[:, 0], xy[:, 1], c=cluster_colors[cluster_labels.cat.codes], s=7, linewidths=0)
    for label in cluster_labels.cat.categories:
        center = np.median(xy[np.asarray(cluster_labels == label)], axis=0)
        ax.text(*center, str(label), fontsize=11, ha='center', va='center',
                bbox={'facecolor':'white', 'edgecolor':'none', 'alpha':.8, 'pad':1})
    ax.set(xlabel='UMAP 1', ylabel='UMAP 2')
    ax.set_title(f'D  Expression | {cluster_labels.nunique()} clusters', loc='left', pad=13)

    ax = fig.add_subplot(grid[1, 1])
    domain_labels = domains.obs.spatial_domain.astype('category')
    domain_colors = plt.get_cmap('tab10')(np.arange(len(domain_labels.cat.categories)))
    polygon_panel(ax, domain_colors[domain_labels.cat.codes])
    handles = [Rectangle((0, 0), 1, 1, facecolor=domain_colors[i], label=str(s))
               for i, s in enumerate(domain_labels.cat.categories)]
    ax.legend(handles=handles, ncol=len(handles), loc='upper center', bbox_to_anchor=(.5, -.12),
              frameon=False, handlelength=1, columnspacing=1, fontsize=10, title='Spatial domain')
    ax.set_title(f'E  Spatial domains | {domain_labels.nunique()} domains', loc='left', pad=13)

    gene = svg.index[0]
    with h5py.File(run/'03_domains/kidney_domains.h5ad', 'r') as handle:
        var = ad.io.read_elem(handle['var'])
        expression = ad.io.read_elem(handle['X'])[:, var.index.get_loc(gene)].toarray().ravel()
    ax = fig.add_subplot(grid[1, 2])
    norm = Normalize(vmin=0, vmax=float(expression.max()))
    polygon_panel(ax, plt.get_cmap('magma')(norm(expression)))
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='magma'), ax=ax, shrink=.8,
                 pad=.025, label='Normalized log1p expression')
    ax.set_title(f'F  Spatial gene | {gene}', loc='left', pad=13)
    ax.text(0, -.14, f"Moran's I = {svg.iloc[0]['I']:.3f}", transform=ax.transAxes, fontsize=11)

    qc_grid = grid[2, 0].subgridspec(1, 2, wspace=.35)
    for i, (column, label) in enumerate([('n_genes_by_counts', 'Genes / cell'), ('transcript_counts', 'Gene transcripts / cell')]):
        ax = fig.add_subplot(qc_grid[0, i])
        values = qc[column].to_numpy()
        positive = values[values > 0]
        ax.hist(positive, bins=np.geomspace(1, max(2, positive.max()), 55), color=['#397AC0', '#E9A31B'][i], edgecolor='white', linewidth=.3)
        median = np.median(positive)
        ax.axvline(median, color='#333333', linestyle='--', linewidth=1.3)
        ax.set_xscale('log'); ax.set(xlabel=label, ylabel='Cells' if i == 0 else '')
        ax.text(.04, .96, f'Median {median:,.0f}', transform=ax.transAxes, va='top', fontsize=11)
        if i == 0:
            ax.set_title('G  Full-sample QC', loc='left', pad=13)
            ax.text(0, -.19, f'{len(qc):,} source cells; distributions show positive counts.', transform=ax.transAxes, fontsize=10, color='#555555')

    ax = fig.add_subplot(grid[2, 1])
    limit = np.nanmax(np.abs(enrichment.to_numpy()))
    heatmap = ax.imshow(enrichment.to_numpy(), cmap='RdBu_r', vmin=-limit, vmax=limit)
    ax.set_xticks(range(len(enrichment.columns)), [abbreviations.get(s, _short_cell_label(s)) for s in enrichment.columns], rotation=55, ha='right')
    ax.set_yticks(range(len(enrichment.index)), [abbreviations.get(s, _short_cell_label(s)) for s in enrichment.index])
    fig.colorbar(heatmap, ax=ax, shrink=.8, pad=.025, label='Z-score')
    ax.set_title('H  Neighborhood enrichment', loc='left', pad=13)
    ax.text(0, -.19, f'{selection["retained_cells"]:,} cells | {len(enrichment)} types', transform=ax.transAxes, fontsize=11)

    ax = fig.add_subplot(grid[2, 2]); ax.set_axis_off()
    weights = communication.groupby(['source', 'target'])['prob'].sum().nlargest(20)
    graph = nx.DiGraph()
    graph.add_nodes_from(enrichment.index)
    graph.add_weighted_edges_from([(source, target, value) for (source, target), value in weights.items()])
    positions = nx.circular_layout(graph)
    max_count = max(counts[s] for s in graph.nodes)
    sizes = [260+900*counts[s]/max_count for s in graph.nodes]
    nx.draw_networkx_nodes(graph, positions, ax=ax, node_color=[palette[s] for s in graph.nodes],
                           node_size=sizes, edgecolors='white', linewidths=1)
    nx.draw_networkx_edges(graph, positions, ax=ax, node_size=sizes, arrows=True, arrowsize=13,
                           width=[.7+4.5*graph[u][v]['weight']/weights.max() for u, v in graph.edges],
                           edge_color=[palette[u] for u, v in graph.edges], alpha=.65,
                           connectionstyle='arc3,rad=0.12')
    for label, position in positions.items():
        ax.text(*(position*1.23), abbreviations.get(label, _short_cell_label(label)),
                ha='center', va='center', fontsize=11)
    ax.set(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), aspect='equal')
    ax.set_title('I  Spatial CellChat network', loc='left', pad=13)
    ax.text(0, -.13, 'Top 20 cell-type pairs by summed LR probability.\nNode size: cells; arrow width: summed score.',
            transform=ax.transAxes, fontsize=10, color='#555555')

    fig.text(.045, .016, 'A-C: official azimuth_coarse annotations | D-F: 2,183-cell ROI | H-I: 1,815 cells after the documented label selection',
             fontsize=12, color='#555555')
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180); plt.close(fig)
    output.with_suffix('.json').write_text(json.dumps({
        'full_cells':len(full), 'roi_cells':loaded.n_obs, 'roi_labels':len(counts),
        'matched_polygons':len(polygons), 'centroid_error_um':error,
        'expression_clusters':int(cluster_labels.nunique()), 'spatial_domains':int(domain_labels.nunique()),
        'svg':str(gene), 'svg_moran_i':float(svg.iloc[0]['I']),
        'downstream_cells':selection['retained_cells'], 'communication_rows':len(communication),
        'displayed_communication_pairs':len(weights), 'figure':str(output),
    }, indent=2), encoding='utf-8')
    return str(output)
