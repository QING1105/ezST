@echo off
REM ezST 一键安装（Windows）— pip 方式
REM 用法: install.bat
REM 1. pip 安装 ezST（自动装全部 Python 依赖）
REM 2. 安装 R 包（Seurat / SPOTlight / CellChat）— 仅 S4/S5 需要
chcp 65001 >nul
setlocal
set SCRIPT_DIR=%~dp0
set REPO_DIR=%SCRIPT_DIR%..

echo === ezST 安装 ===

REM 1. 检查 Python >= 3.10
python -c "import sys; assert sys.version_info >= (3,10)" >nul 2>nul
if errorlevel 1 (
    echo [错误] 需要 Python ^>= 3.10
    pause
    exit /b 1
)

REM 2. pip 安装
echo [1/2] 安装 ezST (pip)...
pip install -e "%REPO_DIR%"
if errorlevel 1 (
    echo [错误] pip 安装失败
    pause
    exit /b 1
)

REM 3. R 包
echo [2/2] 安装 R 包 (Seurat / SPOTlight / CellChat)...
where Rscript >nul 2>nul
if errorlevel 1 (
    echo [警告] 未检测到 R。S4/S5 需要 R ^>= 4.3，请安装后运行: Rscript install\install_r_packages.R
) else (
    Rscript "%SCRIPT_DIR%install_r_packages.R"
)

echo.
echo === ezST 安装完成 ===
echo Python: scanpy/squidpy/scvi-tools/leidenalg 已就绪
echo R: 运行 Rscript install\install_r_packages.R 可重装/更新 R 包
pause
