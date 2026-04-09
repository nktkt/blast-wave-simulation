#!/usr/bin/env pvpython
"""
ParaView 可視化スクリプト
Usage: pvpython paraview_visualize.py
  または ParaView の Python Shell で実行
"""
from paraview.simple import *

# データ読み込み
reader = OpenDataFile("/Users/naoki/dev/lab/love/ver46/results/blast.pvd")
UpdatePipeline()

# カラーマップ設定 (圧力)
display = Show(reader)
ColorBy(display, ('POINTS', 'Pressure'))
pressureLUT = GetColorTransferFunction('Pressure')
pressureLUT.ApplyPreset('Cool to Warm', True)
pressureLUT.RescaleTransferFunction(101325.0, 5066250.0)

# カラーバー
bar = GetScalarBar(pressureLUT)
bar.Title = 'Pressure [Pa]'
bar.Visibility = 1

# 断面スライス (z=0)
slice1 = Slice(Input=reader)
slice1.SliceType = 'Plane'
slice1.SliceType.Origin = [0, 0, 0]
slice1.SliceType.Normal = [0, 0, 1]
slice1Display = Show(slice1)
ColorBy(slice1Display, ('POINTS', 'Pressure'))

# 等値面 (過圧 = 35 kPa, 建物被害閾値)
contour1 = Contour(Input=reader)
contour1.ContourBy = ['POINTS', 'Overpressure']
contour1.Isosurfaces = [7000, 35000, 100000]
contour1Display = Show(contour1)
contour1Display.Opacity = 0.3

# カメラ設定
view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [1920, 1080]
view.Background = [0.1, 0.1, 0.15]
view.CameraPosition = [120, -120, 80]
view.CameraFocalPoint = [0, 0, 10]

# アニメーションコントローラ
scene = GetAnimationScene()
scene.PlayMode = 'Sequence'

Render()
print("ParaView 可視化完了")
print("操作: 再生ボタンでアニメーション, マウスで視点回転")
