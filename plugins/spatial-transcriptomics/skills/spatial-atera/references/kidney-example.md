# 肾脏示例记录

来源：[Atera In Situ Gene Expression, FFPE Human Kidney](https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-nd-kidney)，2026-10-15，CC BY 4.0；实验清单版本 `atera-1.0.4.1`，区域 `S2B_kidney_ND`。

## 数据与运行范围

原始矩阵 399,120 × 24,674，含 18,049 个基因；官方保留 398,873 个细胞。下载脚本取得 `experiment.spatial`、`cell_feature_matrix.zarr.zip`、`cells.zarr.zip` 和 `panel_config.json`，共 1,863,913,974 字节，并核对大小和 SHA256。

局部范围 x=[6300,6900)、y=[3300,3900) µm，保留 2,183 个细胞。分割按 cell ID 对齐，中心坐标差 0 µm。全组织图保留全部 56 类原始标签；区域图包含 24 类，使用同一套颜色。[全组织完整图例](../assets/kidney_full_annotation.png)。

## 复现入口

以下命令在仓库根目录运行，逐阶段查看结果：

```bash
python plugins/spatial-transcriptomics/skills/spatial-atera/scripts/download_kidney.py --output data/atera-kidney
python plugins/spatial-transcriptomics/skills/spatial-atera/scripts/kidney_example.py --stage load --data data/atera-kidney --output results/atera-kidney
```

随后将 `--stage` 依次换成 `expression`、`domains`、`neighborhood`、`communication`。`annotation` 单独重画空间注释；`overview` 从已完成阶段的结果生成展示组图。每阶段输出存于对应的编号子目录。

## 本次参数和结果

| 阶段 | 输入与参数 | 结果 |
|---|---|---|
| 表达 | 区域原始基因计数；总量 10,000、log1p、2,000 HVG、30 PC、Leiden resolution=0.5 | 14 个聚类 |
| 注释 | 官方 `azimuth_coarse`；按 cell ID 对应边界 | 区域 2,183 个细胞、24 类标签 |
| 空间域 / SVG | spatial_weight=0.3、10 个空间近邻、resolution=0.5；Moran's I | 7 个域；首位 CALB1，I≈0.548 |
| 邻域 | 已明确注释且区域内每类至少 10 个细胞；10 个空间近邻 | 1,815 个细胞、14 类，富集 Z-score |
| 通讯 | 同上细胞；CellChat human / spatial；1 µm/坐标单位、tolerance=5、range=250、scale=0.1 | 7,204 条配体—受体记录 |

组图读取保存的分析结果。通讯面板按细胞类型对汇总配体—受体得分，展示得分最高的 20 条有向连线；节点大小表示该类纳入的细胞数。

格式依据：[矩阵与注释](https://www.10xgenomics.com/support/software/atera-onboard-analysis/latest/analysis/output-zarr-cfmatrix)、[分割与坐标](https://www.10xgenomics.com/support/software/atera-onboard-analysis/latest/analysis/output-zarr-cell)。
