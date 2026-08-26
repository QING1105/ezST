#!/bin/bash
# ezST 一键安装（Linux / macOS）— pip 方式
# 用法: bash install.sh
# 1. pip 安装 ezST（自动装全部 Python 依赖）
# 2. 安装 R 包（Seurat / SPOTlight / CellChat）— 仅 S4/S5 需要

set -e
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"

echo -e "${YELLOW}=== ezST 安装 ===${NC}"

# 1. 检查 Python >= 3.10
python3 -c "import sys; assert sys.version_info >= (3,10), '需要 Python >= 3.10'" 2>/dev/null || {
  echo -e "\033[0;31m错误: 需要 Python >= 3.10\033[0m"; exit 1; }

# 2. pip 安装（Python 依赖自动装）
echo -e "${YELLOW}[1/2] 安装 ezST (pip)...${NC}"
pip install -e "$REPO_DIR"

# 3. R 包（仅 S4 反卷积 / S5 通讯需要）
echo -e "${YELLOW}[2/2] 安装 R 包 (Seurat / SPOTlight / CellChat)...${NC}"
if command -v Rscript &> /dev/null; then
  Rscript "$SCRIPT_DIR/install_r_packages.R"
else
  echo -e "\033[0;33m⚠️  未检测到 R。S4/S5 需要 R >= 4.3，请安装后运行: Rscript install/install_r_packages.R\033[0m"
fi

echo -e "\n${GREEN}=== ezST 安装完成 ===${NC}"
echo "Python: scanpy/squidpy/scvi-tools/leidenalg 已就绪"
echo "R: 运行 Rscript install/install_r_packages.R 可重装/更新 R 包"
