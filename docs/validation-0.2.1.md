# 0.2.1 修复与验证

基线：`429f9fc716ed8dde247d72d421d2dfa9ec40cac4`。本地验证日期：2026-10-25。

## 数据与计算方法

自动回归使用固定种子的人工计数与人工坐标。经典流程数据为 120 个位置 × 600 个基因的负二项计数；归一化为总量缩放至 10,000 后取 log1p，同时保存原始计数。Seurat v3 高变基因与 Scanpy 直接调用逐一比较。Moran I 对全部基因计算后排序。原始计数读取顺序为 counts 层、有效 raw.X、整数 X。

空间域测试保持表达不变，只改变坐标；权重 0 的邻接图保持一致，权重 0.5 的邻接图随坐标改变。两张图分别做对称度归一化后混合，再进行 Leiden 聚类。HD 重叠归属测试同时覆盖普通矩阵、稀疏矩阵、等距重叠、未分配 bin 和逐基因计数守恒。

## 本地结果

- 自动回归：31 项通过，0 项失败，约 29.5 秒。执行 `python -m pytest -q`。
- DestVI：CondSCVI 与 DestVI 各训练 1 轮，输出 120 × 2 类型比例；每行比例之和约为 1。
- SPOTlight：实际 R 推理与 CSV/图形导出通过；120 个位置、2 类比例，barcode 正确保留。
- CellChat 2.1.2：90 个细胞 × 15 个配体/受体基因，三类标签。表达模式输出 66 条通讯；空间模式使用相邻坐标时输出 66 条，其中跨类型 40 条。将第二、三类细胞分别沿 x 轴平移 1,000 与 2,000 微米后，距离约束输出 0 条。空间参数为 ratio=1、tol=5 µm、interaction.range=250 µm、scale.distance=1。三种运行均保存 RDS、CSV、运行参数和图形；空结果有列名。
- StarDist：软件包自带 `test_image_he_2d()`，300 × 500 RGB；官方 `2D_versatile_he` 权重；5–95 百分位归一化，prob_thresh=0.01、nms_thresh=0.001、min_area=20。CPU 推理后保留 333 个核，标签与轮廓一致。4 像素网格生成 9,375 个 bin，20 个基因的 Poisson(3) 人工计数，随机种子 25。扩张 5 像素后分配 7,389 个 bin；聚合前后已分配计数均为 443,379，逐基因一致。
- HD 下游：真实运行空间权重 0.5 的域聚类，接 1,000 次标签置换的邻域富集；模块入口、结果读取与导出通过。裁图模块用显式路径完成输出。
- Proseg：依据官方 3.2.0 源码的细胞 × 基因矩阵及元数据格式，验证输出命令参数、基因名称、行列顺序、尺度换算和错误元数据拒绝。此项使用协议测试夹具；完整 Proseg MCMC 计算列入实际样本验收。
- Python 依赖解析通过；R 安装脚本语法检查通过；Python 编译检查、插件运行代码一致性与补丁空白检查通过。

本机环境：Windows、Python 3.11.7、R 4.4.3、Scanpy 1.11.5、Squidpy 1.8.2、scvi-tools 1.4.2、PyTorch 2.14.1+cpu、StarDist 0.9.2、TensorFlow 2.21.0、Seurat 5.4.0、SPOTlight 1.10.0、CellChat 2.1.2。

## 修复范围

依赖约束；原始计数读取与保存；SVG 全基因排名与空间图检查；小鼠线粒体 QC；正确的 Seurat v3 字段；压缩 10x 文件读取；普通矩阵兼容；Windows R 字符串和日志；显式物种；SPOTlight barcode；HD 唯一 bin 归属；分割实例身份；Proseg 基因元数据；可选空间域图；真实 CellChat 空间模式；可移植裁图和下游入口。

补充修复包括：数值 NPZ 几何、Windows StarDist 模型缓存的符号链接权限兼容、零通讯结果导出、下游阶段传递最新数据、scVI 注释原始计数约束、Python/插件版本统一、安装脚本选择当前解释器及验证 R 依赖、插件代码同步工具。

运行 `python -m pytest -q` 执行本地回归。本地记录与合成数据用于软件行为验证，生物学结论在实际样本上结合组织信息、marker 和参考数据审查。

## 方法来源

- [StarDist 官方模型与使用说明](https://github.com/stardist/stardist)
- [CellChat 2.1.2 对象与空间因素](https://github.com/jinworks/CellChat/blob/v2.1.2/R/CellChat_class.R)
- [CellChat 2.1.2 距离计算](https://github.com/jinworks/CellChat/blob/v2.1.2/R/modeling.R)
- [Proseg 官方输出与源码](https://github.com/dcjones/proseg)
