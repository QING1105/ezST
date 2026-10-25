"""Run annotation, domains, neighborhood enrichment and CellChat on HD cells."""
import argparse
import os
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cells', required=True)
    ap.add_argument('--ref', default='')
    ap.add_argument('--cell-type-key', default='cell_type')
    ap.add_argument('--out-dir', default='./hd_downstream')
    ap.add_argument('--n-neighs', type=int, default=6)
    ap.add_argument('--resolution', type=float)
    ap.add_argument('--spatial-weight', type=float, default=0.0)
    for stage in ('annotate', 'domains', 'neighbor', 'cellchat'):
        ap.add_argument(f'--skip-{stage}', action='store_true')
    ap.add_argument('--species', choices=['human', 'mouse'])
    ap.add_argument('--rscript', default='Rscript')
    ap.add_argument('--cellchat-mode', choices=['expression', 'spatial'], default='expression')
    ap.add_argument('--microns-per-coordinate', type=float)
    ap.add_argument('--spatial-tolerance', type=float)
    ap.add_argument('--interaction-range', type=float, default=250.0)
    ap.add_argument('--scale-distance', type=float, default=0.01)
    ap.add_argument('--contact-range', type=float)
    ap.add_argument('--k-min', type=int, default=10)
    args = ap.parse_args()
    import anndata as ad
    try:
        from .visium import identify_spatial_domains, infer_spatial_cell_communication, spatial_neighborhood_enrichment
    except ImportError:
        from visium import identify_spatial_domains, infer_spatial_cell_communication, spatial_neighborhood_enrichment
    os.makedirs(args.out_dir, exist_ok=True)
    current_path = args.cells
    adata = ad.read_h5ad(current_path)
    key = args.cell_type_key
    if 'spatial' not in adata.obsm:
        ap.error("cells requires obsm['spatial']")
    if not args.skip_annotate and key not in adata.obs:
        if not args.ref:
            ap.error('--ref is required to annotate cells')
        script = os.path.join(os.path.dirname(__file__), 'annotate_cells.py')
        subprocess.run([sys.executable, script, '--cells', current_path, '--ref', args.ref,
                        '--cell-type-key', key, '--method', 'auto', '--out-dir', args.out_dir], check=True)
        adata = ad.read_h5ad(os.path.join(args.out_dir, 'annotated_cells.h5ad'))
        if key != 'cell_type':
            adata.obs[key] = adata.obs['cell_type'].copy()
        current_path = os.path.join(args.out_dir, '01_annotated.h5ad')
        adata.write_h5ad(current_path, compression='lzf')
    if not args.skip_domains:
        output = os.path.join(args.out_dir, '02_domains.h5ad')
        print(identify_spatial_domains(current_path, output_path=output,
                                       plot_path=os.path.join(args.out_dir, '02_domains.png'),
                                       resolution=args.resolution, n_neighs=args.n_neighs,
                                       spatial_weight=args.spatial_weight))
        current_path = output
        adata = ad.read_h5ad(current_path)
        print(f"Domains: {adata.obs['spatial_domain'].nunique()}")
    if (not args.skip_neighbor or not args.skip_cellchat) and key not in adata.obs:
        ap.error(f'cells must contain obs[{key!r}] for neighborhood / CellChat')
    if not args.skip_neighbor:
        adata.obs[key] = adata.obs[key].astype('category')
        current_path = os.path.join(args.out_dir, '03_neighborhood_input.h5ad')
        adata.write_h5ad(current_path, compression='lzf')
        print(spatial_neighborhood_enrichment(current_path, cell_type_key=key,
                                             output_path=os.path.join(args.out_dir, '03_neighborhood_enrichment.csv'),
                                             plot_path=os.path.join(args.out_dir, '03_neighborhood.png'),
                                             n_neighs=args.n_neighs))
    if not args.skip_cellchat:
        indices = [group.sample(min(len(group), 1000), random_state=0).index
                   for _, group in adata.obs.groupby(key, observed=True)]
        import numpy as np
        cc = adata[np.concatenate(indices)].copy()
        cc_path = os.path.join(args.out_dir, '04_cellchat_input.h5ad')
        cc.write_h5ad(cc_path, compression='lzf')
        print(infer_spatial_cell_communication(
            cc_path, cell_type_key=key, output_dir=os.path.join(args.out_dir, 'cellchat_results'),
            species=args.species, r_script_path=args.rscript, mode=args.cellchat_mode,
            microns_per_coordinate=args.microns_per_coordinate, spatial_tolerance=args.spatial_tolerance,
            interaction_range=args.interaction_range, scale_distance=args.scale_distance,
            contact_range=args.contact_range, k_min=args.k_min))
    print('DONE')


if __name__ == '__main__':
    main()
