# 把 patent-lens-to-excel 推到 GitHub（在本机 PowerShell 里跑，用本机 git 凭证 / Windows 凭据管理器）
# 用法：powershell -ExecutionPolicy Bypass -File tools\push-to-github.ps1 "提交说明"
# 仓库：https://github.com/anvcor/patent-lens-to-excel（main 分支，直推）
param([string]$Message = "更新 patent-lens-to-excel")
$ErrorActionPreference = 'Stop'
$src = Split-Path -Parent $PSScriptRoot          # tools\ 的上一级 = 仓库根 E:\Download\patent-lens-to-excel
Set-Location $src
git add -A
git status --short
# 全局没配 user.name/email，提交时临时给（与 cv2zmx、GlassMapTool 两个仓库一致）
git -c user.name="T1791" -c user.email="t1791669914@gmail.com" commit -m $Message
git push origin main
