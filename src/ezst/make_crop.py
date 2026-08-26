#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""make_crop.py — 用 ROI 坐标生成高清 crop 图

用法:
  python make_crop.py --method stardist|proseg --roi x0,y0,x1,y1 [--out out.png]
  python make_crop.py --method proseg --roi-json crops/roi_<sample>_proseg.json

输出统一到 results/07b_segmentation/crops/ 子目录。
"""
import argparse, json, os, sys
sys.path.insert(0, '/root/scripts/visium_hd')
import visium_new_platforms as vnp
import anndata as ad

HE = '/root/hest_data/visium_hd_official/tissue_image.btf'
BASE = '/root/scripts/visium_hd/results/07b_segmentation'
CROPS_DIR = os.path.join(BASE, 'crops')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--method', choices=['stardist', 'proseg'], required=True)
    ap.add_argument('--sample', default='visium_hd_official',
                    help='样本名 (用于输出文件命名)')
    ap.add_argument('--roi', help='x0,y0,x1,y1 (全分辨率像素)')
    ap.add_argument('--roi-json', help='pick_roi.py 输出的 json 文件')
    ap.add_argument('--out', default=None, help='自定义输出路径 (默认 crops/)')
    ap.add_argument('--size', type=int, default=None,
                    help='若 roi 是中心点, 用 size 扩展')
    args = ap.parse_args()

    if args.roi_json:
        with open(args.roi_json) as f:
            r = json.load(f)
        x0, y0, x1, y1 = r['x0'], r['y0'], r['x1'], r['y1']
    elif args.roi:
        parts = [int(v) for v in args.roi.split(',')]
        if len(parts) == 4:
            x0, y0, x1, y1 = parts
        elif len(parts) == 2 and args.size:
            cx, cy = parts
            s = args.size // 2
            x0, y0, x1, y1 = cx - s, cy - s, cx + s, cy + s
        else:
            raise SystemExit('--roi 需要 4 个值, 或 2 个值+--size')
    else:
        raise SystemExit('需要 --roi 或 --roi-json')

    if args.method == 'stardist':
        a = ad.read_h5ad(f'{BASE}/{args.sample}_stardist/04_bin2cell.h5ad')
        geom = f'{BASE}/{args.sample}_stardist/03_nuclei_polys.pckl'
        color = 'tab:blue'
        kind = 'stardist'
    else:
        a = ad.read_h5ad(f'{BASE}/{args.sample}_proseg/04_proseg_cells.h5ad')
        geom = f'{BASE}/{args.sample}_proseg/proseg_out.zarr/shapes/cell_boundaries/shapes.parquet'
        color = 'tab:red'
        kind = 'proseg'

    # 默认输出到 crops/ 子目录 (自动创建)
    os.makedirs(CROPS_DIR, exist_ok=True)
    out = args.out or os.path.join(
        CROPS_DIR,
        f'crop_{args.sample}_{args.method}_{x0}_{y0}_{x1}_{y1}.png')
    vnp.plot_cells_crop(
        a, HE, out,
        roi=(x0, y0, x1, y1),
        geometry_path=geom, geometry_kind=kind,
        edge_color=color, linewidth=1.2, dpi=150,
        title=f'{args.method} {args.sample} crop x{x0}:{x1} y{y0}:{y1}')
    print(f'SAVED {out}', flush=True)


if __name__ == '__main__':
    main()
