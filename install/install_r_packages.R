# install_r_packages.R — ezST R 依赖安装
# 用法: Rscript install_r_packages.R
# 安装: Seurat 系列 + SPOTlight (Bioconductor) + CellChat v2 (GitHub)

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)  # 大包下载超时放宽

cat("=== ezST R 依赖安装 ===\n")

# 0. 检查编译工具（Windows 需要 Rtools）
if (.Platform$OS.type == "windows") {
  rt <- Sys.getenv("RTOOLS_VERSION")
  if (rt == "") {
    cat("⚠️  未检测到 Rtools。如果安装需要编译的包失败，请先安装 Rtools:\n")
    cat("    https://cran.r-project.org/bin/windows/Rtools/\n")
  } else {
    cat("✓ 检测到 Rtools", rt, "\n")
  }
} else {
  cat("✓ Linux/macOS: 编译工具由系统包管理器提供\n")
}

# 1. 基础工具
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# 2. Seurat 系列（S5 CellChat 需要 Seurat 对象）
for (pkg in c("Seurat", "SeuratObject")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("安装", pkg, "...\n")
    install.packages(pkg)
  } else {
    cat("✓", pkg, "已安装\n")
  }
}

# 3. SPOTlight（S4 反卷积，Bioconductor）
if (!requireNamespace("SPOTlight", quietly = TRUE)) {
  cat("安装 SPOTlight (Bioconductor)...\n")
  BiocManager::install("SPOTlight", update = FALSE, ask = FALSE)
} else {
  cat("✓ SPOTlight 已安装\n")
}

# SPOTlight 依赖的 Bioconductor 基础包
for (pkg in c("SingleCellExperiment", "SummarizedExperiment", "scater", "scran")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("安装", pkg, "(Bioconductor)...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# 4. CellChat v2（S5 细胞通讯，仅 GitHub）
if (!requireNamespace("CellChat", quietly = TRUE)) {
  cat("安装 CellChat v2 (GitHub jinworks/CellChat)...\n")
  devtools::install_github("jinworks/CellChat")
} else {
  cat("✓ CellChat 已安装\n")
}

cat("\n=== ezST R 依赖安装完成 ===\n")
cat("验证: Rscript -e 'library(Seurat); library(SPOTlight); library(CellChat)'\n")
