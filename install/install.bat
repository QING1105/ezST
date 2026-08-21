@echo off
REM ezST 一键环境安装（Windows）
REM 用法: install.bat
REM 创建 conda 环境 ezST → 安装 Python 包 → 安装 R 包 → 验证
chcp 65001 >nul
setlocal enabledelayedexpansion

set ENV_NAME=ezST
set SCRIPT_DIR=%~dp0

echo === ezST 环境安装 ===

REM 1. 检查 conda
where conda >nul 2>nul
if errorlevel 1 (
    echo [错误] 未检测到 conda。请先安装 Miniconda: https://docs.conda.io/en/latest/miniconda.html
    pause
    exit /b 1
)

REM 2. 创建/更新 conda 环境
conda env list | findstr /r /c:"^ezST " >nul
if errorlevel 1 (
    echo [创建] 环境 %ENV_NAME% ...
    conda env create -n %ENV_NAME% -f "%SCRIPT_DIR%environment.yml"
) else (
    echo [更新] 环境 %ENV_NAME% ...
    conda env update -n %ENV_NAME% -f "%SCRIPT_DIR%environment.yml"
)
if errorlevel 1 (
    echo [错误] conda 环境创建失败
    pause
    exit /b 1
)

REM 3. 激活环境
call conda activate %ENV_NAME%
echo [✓] 环境 %ENV_NAME% 已激活

REM 4. pip 补充依赖
echo [安装] Python 补充依赖...
pip install -r "%SCRIPT_DIR%requirements.txt"

REM 5. R 包安装
echo [安装] R 包 (Seurat / SPOTlight / CellChat)...
Rscript "%SCRIPT_DIR%install_r_packages.R"

REM 6. 验证
echo [验证] 安装结果...
python -c "import scanpy, squidpy, scvi; print('[OK] scanpy', scanpy.__version__); print('[OK] squidpy', squidpy.__version__); print('[OK] scvi-tools', scvi.__version__)"
Rscript -e "library(Seurat); library(SPOTlight); library(CellChat); cat('[OK] Seurat / SPOTlight / CellChat\n')"

echo.
echo === ezST 环境安装完成 ===
echo 使用: conda activate %ENV_NAME%
pause
