# -*- coding: utf-8 -*-
"""
pick_roi.py — Visium HD 交互式 ROI 选择工具

用法:
    python pick_roi.py --image tissue_image.png [--scale 0.2738] [--out roi.json]

功能:
    1. 打开图片 (自动降采样到适合屏幕, 记录缩放因子)
    2. 鼠标左键拖拽框选 ROI (松开确认, 可多次重选)
    3. 键盘: [Enter]=确认输出  [r]=重选  [q]=退出
    4. 输出全分辨率像素坐标到 json: {"x0":..,"y0":..,"x1":..,"y1":..,"scale":..}
    5. 若提供 --save-crop, 保存一张 ROI 局部预览图

依赖: pip install matplotlib pillow numpy
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector


def load_image(path, max_display=1600):
    """读图并降采样到适合屏幕; 返回 (img, scale_factor 显示->原始)."""
    from PIL import Image
    img = Image.open(path)
    if img.mode == "RGBA":
        img = img.convert("RGB")
    w, h = img.size
    scale = 1.0
    if max(w, h) > max_display:
        scale = max_display / max(w, h)
        img = img.resize((int(w * scale), int(h * scale)), Image.LANCZOS)
    return np.asarray(img), scale


class ROISelector:
    def __init__(self, img, scale, out_path, save_crop=None):
        self.img = img
        self.scale = scale  # 显示尺寸 -> 原始尺寸 的乘数 (显示*scale=原始)
        self.out_path = out_path
        self.save_crop = save_crop
        self.roi_disp = None  # 显示坐标 (x0,y0,x1,y1)

        self.fig, self.ax = plt.subplots(figsize=(12, 10))
        self.ax.imshow(img)
        self.ax.set_title("拖拽鼠标框选 ROI | Enter=确认  r=重选  q=退出")
        self.rs = RectangleSelector(
            self.ax, self.on_select,
            interactive=True, useblit=True,
            props=dict(facecolor="cyan", edgecolor="red", alpha=0.25, fill=True),
        )
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        plt.show()

    def on_select(self, eclick, erelease):
        x0, y0 = eclick.xdata, eclick.ydata
        x1, y1 = erelease.xdata, erelease.ydata
        self.roi_disp = (int(min(x0, x1)), int(min(y0, y1)),
                         int(max(x0, x1)), int(max(y0, y1)))
        print(f"[ROI 显示坐标] x[{self.roi_disp[0]}:{self.roi_disp[2]}] "
              f"y[{self.roi_disp[1]}:{self.roi_disp[3]}]", flush=True)

    def on_key(self, event):
        if event.key == "enter":
            self.confirm()
        elif event.key == "r":
            self.roi_disp = None
            self.ax.set_title("已重选: 重新拖拽框选")
            self.fig.canvas.draw_idle()
        elif event.key == "q":
            print("已退出 (未保存)", flush=True)
            plt.close(self.fig)

    def confirm(self):
        if self.roi_disp is None:
            print("还没框选区域, 请先拖拽选择", flush=True)
            return
        x0, y0, x1, y1 = self.roi_disp
        # 显示坐标 -> 全分辨率原始坐标
        ox0, oy0 = int(x0 / self.scale), int(y0 / self.scale)
        ox1, oy1 = int(x1 / self.scale), int(y1 / self.scale)
        roi = {"x0": ox0, "y0": oy0, "x1": ox1, "y1": oy1,
               "display_scale": self.scale}
        with open(self.out_path, "w", encoding="utf-8") as f:
            json.dump(roi, f, indent=2)
        print(f"[已保存] {self.out_path}", flush=True)
        print(f"[全分辨率 ROI] x[{ox0}:{ox1}] y[{oy0}:{oy1}]", flush=True)
        if self.save_crop:
            crop = self.img[y0:y1, x0:x1]
            from PIL import Image
            Image.fromarray(crop).save(self.save_crop)
            print(f"[crop 预览] {self.save_crop}", flush=True)
        plt.close(self.fig)


def main():
    ap = argparse.ArgumentParser(description="Visium HD 交互式 ROI 选择")
    ap.add_argument("--image", required=True, help="HE 图片路径 (png/jpg/tif)")
    ap.add_argument("--scale", type=float, default=None,
                    help="微米/像素 (仅用于展示, 非必需)")
    ap.add_argument("--out", default="roi.json", help="输出坐标 json 路径")
    ap.add_argument("--save-crop", default=None, help="保存 ROI 预览图路径")
    args = ap.parse_args()

    img, scale = load_image(args.image)
    print(f"图片加载完成: 显示 {img.shape[1]}x{img.shape[0]}, "
          f"缩放因子={scale:.3f}", flush=True)
    ROISelector(img, scale, args.out, args.save_crop)


if __name__ == "__main__":
    main()
