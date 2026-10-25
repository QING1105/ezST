"""Create a cell-boundary crop from explicitly supplied sample files."""
import argparse
import json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--method', choices=['stardist', 'proseg'], required=True)
    ap.add_argument('--image', required=True, help='Full-resolution tissue image')
    ap.add_argument('--cells', required=True, help='Cell-level h5ad')
    ap.add_argument('--geometry', required=True, help='StarDist NPZ or Proseg shapes.parquet')
    ap.add_argument('--sample', default='sample')
    roi = ap.add_mutually_exclusive_group(required=True)
    roi.add_argument('--roi', help='x0,y0,x1,y1 or cx,cy with --size (full-resolution pixels)')
    roi.add_argument('--roi-json', help='ROI JSON produced by pick_roi.py')
    ap.add_argument('--out', default=None)
    ap.add_argument('--size', type=int)
    ap.add_argument('--allow-legacy-pickle', action='store_true', help='Load trusted legacy StarDist pickle')
    args = ap.parse_args()
    if args.roi_json:
        with open(args.roi_json, encoding='utf-8') as handle:
            data = json.load(handle)
        bounds = tuple(int(data[key]) for key in ('x0', 'y0', 'x1', 'y1'))
    else:
        parts = [int(value) for value in args.roi.split(',')]
        if len(parts) == 4:
            bounds = tuple(parts)
        elif len(parts) == 2 and args.size and args.size > 0:
            cx, cy = parts
            half = args.size // 2
            bounds = (cx - half, cy - half, cx - half + args.size, cy - half + args.size)
        else:
            ap.error('--roi requires four values or two values with a positive --size')
    x0, y0, x1, y1 = bounds
    if x1 <= x0 or y1 <= y0:
        ap.error('ROI must have positive width and height')
    import anndata as ad
    try:
        from .visium_new_platforms import plot_cells_crop
    except ImportError:
        from visium_new_platforms import plot_cells_crop
    output = Path(args.out) if args.out else Path(args.cells).parent / 'crops' / f'crop_{args.sample}_{args.method}_{x0}_{y0}_{x1}_{y1}.png'
    output.parent.mkdir(parents=True, exist_ok=True)
    plot_cells_crop(ad.read_h5ad(args.cells), args.image, str(output), roi=bounds,
                    geometry_path=args.geometry, geometry_kind=args.method,
                    allow_legacy_pickle=args.allow_legacy_pickle,
                    edge_color='tab:blue' if args.method == 'stardist' else 'tab:red',
                    linewidth=1.2, dpi=150, title=f'{args.method} {args.sample}')
    print(f'SAVED {output}', flush=True)


if __name__ == '__main__':
    main()
