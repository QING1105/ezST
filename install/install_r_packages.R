# install_r_packages.R — ezST R 依赖安装
# 用法: Rscript install_r_packages.R
# 安装: Seurat 系列 + SPOTlight (Bioconductor) + CellChat v2.1.2 (GitHub 稳定版)
#
# 2026-08-26 修订 (踩坑修复):
#   1. CellChat 改装 GitHub release v2.1.2 (稳定版), 不再装 main 开发版 (2.2.0.9001 lazy-load bug)
#   2. 补全 CellChat 依赖链 (ggpubr/doBy/pbkrtest/car/rstatix), 之前缺失导致安装失败
#   3. conda 版 R 需把 <env>/bin 加入 PATH, 否则找不到 x86_64-conda-linux-gnu-cc 编译器
#   4. CellChat 源码包自带预编译 .o/.so (作者环境), 必须删除强制重编, 否则 "无效的 ELF 头"
#   5. v2.1.2 lazy-load 正常: 直接 CellChatDB.mouse 即可, 不要加 data() 显式加载 (会破坏 lazy-load)

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)  # 大包下载超时放宽

# 0. conda 版 R 编译器 PATH 修复 (Linux/macOS)
if (.Platform$OS.type != "windows") {
  # 如果 R 是 conda 装的, 把 conda env bin 加到 PATH 最前 (编译器 x86_64-conda-linux-gnu-cc)
  conda_env <- Sys.getenv("CONDA_PREFIX")
  if (nzchar(conda_env) && dir.exists(file.path(conda_env, "bin"))) {
    Sys.setenv(PATH = paste(file.path(conda_env, "bin"), Sys.getenv("PATH"), sep = ":"))
    cat("✓ conda R 检测到, PATH 已加", file.path(conda_env, "bin"), "\n")
  }
}

cat("=== ezST R 依赖安装 ===\n")

# 1. 检查编译工具
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
  # 验证 conda 编译器可见
  if (nzchar(Sys.getenv("CONDA_PREFIX"))) {
    cc <- Sys.which("x86_64-conda-linux-gnu-cc")
    if (nzchar(cc)) cat("✓ conda 编译器:", cc, "\n")
    else cat("⚠️  x86_64-conda-linux-gnu-cc 未找到, 请确认 conda env 里有编译器 (conda install -c conda-forge cxx-compiler)\n")
  }
}

# 2. 基础工具
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# 3. Seurat 系列（CellChat 需要 Seurat 对象）
for (pkg in c("Seurat", "SeuratObject")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("安装", pkg, "...\n")
    install.packages(pkg)
  } else {
    cat("✓", pkg, "已安装\n")
  }
}

# 4. SPOTlight（反卷积，Bioconductor）
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

# 5. CellChat v2.1.2 稳定版 (S5 细胞通讯)
#    GitHub release: https://github.com/jinworks/CellChat/archive/refs/tags/v2.1.2.tar.gz
#    注意: 不装 main 开发版 (2.2.0.9001), 它有 lazy-load bug 导致 CellChatDB 为空
CELLCHAT_VER <- "2.1.2"
if (!requireNamespace("CellChat", quietly = TRUE)) {
  cat("安装 CellChat v", CELLCHAT_VER, " (GitHub 稳定版)...\n", sep = "")
  # 5.1 先装依赖 (CellChat 的 NAMESPACE 硬引用 ggpubr, 必须先装)
  deps <- c("Deriv", "doBy", "pbkrtest", "car", "rstatix", "ggpubr")
  for (p in deps) {
    if (!requireNamespace(p, quietly = TRUE)) {
      cat("  安装依赖", p, "...\n")
      install.packages(p)
    }
  }
  # 5.2 下载 v2.1.2 源码
  tmp <- tempdir()
  tarball <- file.path(tmp, paste0("CellChat-", CELLCHAT_VER, ".tar.gz"))
  url <- sprintf("https://github.com/jinworks/CellChat/archive/refs/tags/v%s.tar.gz", CELLCHAT_VER)
  download.file(url, tarball, quiet = TRUE)
  untar(tarball, exdir = tmp)
  srcdir <- file.path(tmp, sprintf("CellChat-%s", CELLCHAT_VER))
  # 5.3 删除源码里预编译的 .o/.so (作者环境编译的, 直接装会 "无效的 ELF 头")
  obj_files <- list.files(file.path(srcdir, "src"), pattern = "\\.(o|so)$", full.names = TRUE)
  if (length(obj_files) > 0) {
    file.remove(obj_files)
    cat("  ✓ 已删除", length(obj_files), "个预编译 .o/.so, 强制重新编译\n")
  }
  # 5.4 编译安装
  install.packages(srcdir, repos = NULL, type = "source")
} else {
  # 已装但可能是开发版 -> 提示
  ver <- as.character(utils::packageVersion("CellChat"))
  cat("✓ CellChat 已安装 (版本", ver, ")\n")
  if (ver != CELLCHAT_VER) {
    cat("⚠️  当前版本", ver, "不是稳定版", CELLCHAT_VER, "。如需稳定版请先 remove.packages('CellChat') 再重跑本脚本\n")
  }
}

cat("\n=== ezST R 依赖安装完成 ===\n")
cat("验证: Rscript -e 'library(Seurat); library(SPOTlight); library(CellChat); data(\"CellChatDB.mouse\"); dim(CellChatDB.mouse$interaction)'\n")
