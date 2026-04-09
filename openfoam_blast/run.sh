#!/bin/bash
# OpenFOAM 爆風シミュレーション実行スクリプト
# Docker 内で実行される

set -e

echo "=========================================="
echo "OpenFOAM 爆風シミュレーション (建物あり)"
echo "=========================================="

# 1. メッシュ生成
echo "[1/4] blockMesh - ベースメッシュ生成..."
blockMesh

# 2. 建物メッシュ埋め込み
echo "[2/4] snappyHexMesh - 建物をメッシュに追加..."
snappyHexMesh -overwrite

# 3. 初期条件設定 (爆発源)
echo "[3/4] setFields - 爆発源の高圧球を設定..."
setFields

# 4. ソルバー実行
echo "[4/4] rhoCentralFoam - 圧縮性流れ計算..."
echo "  ドメイン: 200m x 200m x 80m"
echo "  建物: 4棟"
echo "  終了時刻: 0.5 s"
rhoCentralFoam

echo "=========================================="
echo "計算完了！"
echo "ParaView で結果を確認:"
echo "  paraview openfoam_blast/case.foam"
echo "=========================================="
