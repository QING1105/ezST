#!/bin/bash
# ezST 一键环境安装（Linux / macOS）
# 用法: bash install.sh
# 创建 conda 环境 ezST → 安装 Python 包 → 安装 R 包 → 验证

set -e

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

ENV_NAME="ezST"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo -e "${YELLOW}=== ezST 环境安装 ===${NC}"

# 1. 检查 conda
if ! command -v conda &> /dev/null; then
  echo -e "${RED}错误: 未检测到 conda。请先安装 Miniconda: https://docs.conda.io/en/latest/miniconda.html${NC}"
  exit 1
fi

# 2. 创建/更新 conda 环境
if conda env list | grep -q "^${ENV_NAME} "; then
  echo -e "${YELLOW}环境 ${ENV_NAME} 已存在，更新...${NC}"
  conda env update -n "${ENV_NAME}" -f "${SCRIPT_DIR}/environment.yml"
else
  echo -e "${YELLOW}创建环境 ${ENV_NAME}...${NC}"
  conda env create -n "${ENV_NAME}" -f "${SCRIPT_DIR}/environment.yml"
fi

# 3. 激活环境
eval "$(conda shell.bash hook)"
conda activate "${ENV_NAME}"
echo -e "${GREEN}✓ 环境 ${ENV_NAME} 已激活${NC}"

# 4. pip 补充依赖
echo -e "${YELLOW}安装 Python 补充依赖...${NC}"
pip install -r "${SCRIPT_DIR}/requirements.txt"

# 5. R 包安装
echo -e "${YELLOW}安装 R 包 (Seurat / SPOTlight / CellChat)...${NC}"
Rscript "${SCRIPT_DIR}/install_r_packages.R"

# 6. 验证
echo -e "${YELLOW}验证安装...${NC}"
python -c "import scanpy, squidpy, scvi; print('✓ scanpy', scanpy.__version__); print('✓ squidpy', squidpy.__version__); print('✓ scvi-tools', scvi.__version__)"
Rscript -e "library(Seurat); library(SPOTlight); library(CellChat); cat('✓ Seurat / SPOTlight / CellChat OK\n')"

echo -e "\n${GREEN}=== ezST 环境安装完成 ===${NC}"
echo -e "使用: conda activate ${ENV_NAME}"
