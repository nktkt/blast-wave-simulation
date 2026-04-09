#!/usr/bin/env pvpython
"""
ParaView で OpenFOAM 爆風シミュレーションのアニメーション GIF を作成

Usage:
  /Applications/ParaView-6.1.0.app/Contents/bin/pvpython make_gif.py
"""

from paraview.simple import *
import os

# =============================================================================
# 設定
# =============================================================================
CASE_FOAM = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'openfoam_blast', 'case.foam')
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
FRAME_DIR = os.path.join(OUTPUT_DIR, 'gif_frames')
GIF_PATH = os.path.join(OUTPUT_DIR, 'blast_openfoam.gif')
WIDTH = 960
HEIGHT = 720

os.makedirs(FRAME_DIR, exist_ok=True)

# =============================================================================
# データ読み込み
# =============================================================================
print("データ読み込み中...")
reader = OpenFOAMReader(FileName=CASE_FOAM)
reader.MeshRegions = ['internalMesh']
reader.CellArrays = ['p', 'T', 'U']
UpdatePipeline()

# タイムステップ取得
timesteps = reader.TimestepValues
print(f"タイムステップ: {len(timesteps)} 個")
print(f"  t = {timesteps[0]:.4f} ~ {timesteps[-1]:.4f} s")

# =============================================================================
# 断面スライス (z = 5m)
# =============================================================================
print("スライス作成中...")
slice1 = Slice(Input=reader)
slice1.SliceType = 'Plane'
slice1.SliceType.Origin = [0.0, 0.0, 5.0]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
UpdatePipeline()

# =============================================================================
# 表示設定
# =============================================================================
view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [WIDTH, HEIGHT]
view.Background = [0.1, 0.1, 0.15]

# スライスの表示
display = Show(slice1, view)
ColorBy(display, ('CELLS', 'p'))
display.SetRepresentationType('Surface')

# カラーマップ設定
pLUT = GetColorTransferFunction('p')
pLUT.ApplyPreset('Cool to Warm (Extended)', True)
pLUT.RescaleTransferFunction(80000, 500000)

# カラーバー
colorBar = GetScalarBar(pLUT, view)
colorBar.Title = 'Pressure [Pa]'
colorBar.ComponentTitle = ''
colorBar.Visibility = 1
colorBar.TitleFontSize = 14
colorBar.LabelFontSize = 12

# 元データも輪郭だけ表示
outlineDisplay = Show(reader, view)
outlineDisplay.SetRepresentationType('Outline')
outlineDisplay.AmbientColor = [0.5, 0.5, 0.5]
outlineDisplay.DiffuseColor = [0.5, 0.5, 0.5]

# カメラ設定 (上から見下ろし + 少し角度)
view.CameraPosition = [0, -180, 150]
view.CameraFocalPoint = [0, 0, 5]
view.CameraViewUp = [0, 0, 1]
view.CameraParallelScale = 120

Render()

# =============================================================================
# フレーム書き出し
# =============================================================================
print(f"フレーム書き出し中 ({len(timesteps)} 枚)...")

scene = GetAnimationScene()
scene.PlayMode = 'Snap To TimeSteps'

for i, t in enumerate(timesteps):
    scene.AnimationTime = t
    UpdatePipeline(time=t)

    # タイトルにタイムスタンプを表示
    view.OrientationAxesVisibility = 0

    Render()

    frame_path = os.path.join(FRAME_DIR, f'frame_{i:04d}.png')
    SaveScreenshot(frame_path, view,
                   ImageResolution=[WIDTH, HEIGHT],
                   TransparentBackground=0)
    print(f"  [{i+1}/{len(timesteps)}] t={t:.4f}s → {os.path.basename(frame_path)}")

print(f"\n全 {len(timesteps)} フレーム書き出し完了")
print(f"フレーム保存先: {FRAME_DIR}/")

# =============================================================================
# GIF 作成 (Pillow)
# =============================================================================
print("\nGIF作成中...")
try:
    from PIL import Image, ImageDraw, ImageFont

    frames = []
    frame_files = sorted([f for f in os.listdir(FRAME_DIR) if f.endswith('.png')])

    for i, fname in enumerate(frame_files):
        img = Image.open(os.path.join(FRAME_DIR, fname))

        # タイムスタンプを画像に描画
        draw = ImageDraw.Draw(img)
        t = timesteps[i] if i < len(timesteps) else 0
        text = f"t = {t*1000:.1f} ms"
        draw.text((20, 20), text, fill=(255, 255, 255))

        frames.append(img)

    # GIF保存 (各フレーム200ms, ループ)
    frames[0].save(
        GIF_PATH,
        save_all=True,
        append_images=frames[1:],
        duration=200,
        loop=0,
        optimize=True,
    )
    gif_size = os.path.getsize(GIF_PATH) / 1024 / 1024
    print(f"GIF保存: {GIF_PATH} ({gif_size:.1f} MB)")

except ImportError:
    print("Pillow未インストール — ffmpegで作成します")
    import subprocess
    subprocess.run([
        'ffmpeg', '-y', '-framerate', '5',
        '-i', os.path.join(FRAME_DIR, 'frame_%04d.png'),
        '-vf', 'palettegen',
        os.path.join(FRAME_DIR, 'palette.png')
    ], capture_output=True)
    subprocess.run([
        'ffmpeg', '-y', '-framerate', '5',
        '-i', os.path.join(FRAME_DIR, 'frame_%04d.png'),
        '-i', os.path.join(FRAME_DIR, 'palette.png'),
        '-lavfi', 'paletteuse',
        GIF_PATH
    ], capture_output=True)
    gif_size = os.path.getsize(GIF_PATH) / 1024 / 1024
    print(f"GIF保存: {GIF_PATH} ({gif_size:.1f} MB)")

print("\n完了！")
