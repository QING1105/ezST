"""Run one reviewable stage of the official Atera human kidney example."""
import argparse
from pathlib import Path
import json
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', required=True, choices=['load', 'expression', 'annotation', 'domains', 'neighborhood', 'communication', 'overview'])
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--roi', type=float, nargs=4, default=[6300, 3300, 6900, 3900],
                        metavar=('XMIN', 'YMIN', 'XMAX', 'YMAX'), help='Bounds in micrometres')
    parser.add_argument('--full', action='store_true', help='Load all officially retained cells instead of the demo ROI')
    parser.add_argument('--rscript', default='Rscript')
    args = parser.parse_args()
    try:
        from ezst.atera import load_atera_zarr, plot_atera_qc, plot_atera_segmentation, plot_atera_annotations
        from ezst import visium
    except ImportError:
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'spatial-visium-hd'))
        from atera import load_atera_zarr, plot_atera_qc, plot_atera_segmentation, plot_atera_annotations
        import visium
    out, data = args.output.resolve(), args.data.resolve()
    out.mkdir(parents=True, exist_ok=True)
    loaded = out/'01_loading/kidney_cells.h5ad'
    normalized = out/'02_expression/kidney_normalized.h5ad'
    clustered = out/'02_expression/kidney_clustered.h5ad'
    domains = out/'03_domains/kidney_domains.h5ad'
    stages = {'load':'01_loading', 'expression':'02_expression', 'annotation':'02_expression', 'domains':'03_domains',
              'neighborhood':'04_neighborhood', 'communication':'05_communication', 'overview':'06_overview'}
    current = out/stages[args.stage]
    current.mkdir(exist_ok=True)
    logs = []
    if args.stage == 'load':
        logs.append(load_atera_zarr(data, loaded, sample_id='Atera_human_kidney',
                                    roi=None if args.full else args.roi))
        logs.append(plot_atera_qc(loaded, current))
        if not args.full:
            logs.append(json.dumps(plot_atera_segmentation(loaded, data/'cells.zarr.zip', current/'segmentation.png')))
    elif args.stage == 'expression':
        logs.append(visium.normalize_visium(str(loaded), str(normalized), n_top_genes=2000))
        logs.append(visium.cluster_spatial_data(str(normalized), str(clustered), n_pcs=30, resolution=.5))
        logs.append(plot_atera_annotations(loaded, current/'official_annotations.png',
                     cells_path=data/'cells.zarr.zip', source_matrix_path=data/'cell_feature_matrix.zarr.zip', title='Atera human kidney'))
    elif args.stage == 'annotation':
        logs.append(plot_atera_annotations(loaded, current/'official_annotations.png',
                     cells_path=data/'cells.zarr.zip', source_matrix_path=data/'cell_feature_matrix.zarr.zip', title='Atera human kidney'))
    elif args.stage == 'overview':
        from render_overview import render_overview
        logs.append(render_overview(data, out, current/'kidney_overview.png'))
    elif args.stage == 'domains':
        logs.append(visium.identify_spatial_domains(str(normalized), str(domains), resolution=.5,
                                                    spatial_weight=.3, n_neighs=10, n_pcs=30, n_top_genes=2000))
        logs.append(visium.find_spatially_variable_genes(str(domains), str(current/'svg.csv'), n_top_genes=2000))
    else:
        import anndata as ad
        a = ad.read_h5ad(loaded)
        labels = a.obs['cell_type'].astype(str)
        counts = labels.value_counts()
        keep = ~labels.str.contains('Not Confidently Annotated|Unassigned|Unannotated', regex=True)
        keep &= labels.map(counts).ge(10)
        report = {'input_cells': a.n_obs, 'retained_cells': int(keep.sum()),
                  'excluded_cells': int((~keep).sum()), 'label_source': 'azimuth_coarse',
                  'rule': 'Exclude unresolved official labels and groups with fewer than 10 cells in the ROI'}
        a = a[keep].copy()
        a.obs['cell_type'] = a.obs['cell_type'].cat.remove_unused_categories()
        downstream = current/'analysis_cells.h5ad'
        a.write_h5ad(downstream, compression='lzf')
        (current/'cell_selection.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        del a
        if args.stage == 'neighborhood':
            logs.append(visium.spatial_neighborhood_enrichment(str(downstream), 'cell_type',
                                                              str(current/'enrichment.csv'), n_neighs=10))
        else:
            logs.append(visium.infer_spatial_cell_communication(str(downstream), 'cell_type', str(current),
                          r_script_path=args.rscript, species='human', mode='spatial', microns_per_coordinate=1,
                          spatial_tolerance=5, interaction_range=250, scale_distance=.1))
    text = '\n'.join(logs)
    log_name = 'annotation.log' if args.stage == 'annotation' else 'run.log'
    (current/log_name).write_text(text, encoding='utf-8')
    print(text, flush=True)


if __name__ == '__main__':
    main()
