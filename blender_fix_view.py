"""
Blender で blast_wave.blend を開いた状態で実行するビュー修正スクリプト
Scripting タブで実行、またはコマンドラインで:
  blender results/blast_wave.blend --python blender_fix_view.py
"""
import bpy
import os

# 全ての 3D ビューポートをマテリアルプレビューに変更
for area in bpy.context.screen.areas:
    if area.type == 'VIEW_3D':
        for space in area.spaces:
            if space.type == 'VIEW_3D':
                space.shading.type = 'MATERIAL'  # マテリアルプレビュー
                space.clip_end = 1000

# フレーム5に設定 (爆発が見える時点)
bpy.context.scene.frame_set(5)

# カメラビューに切り替え
for area in bpy.context.screen.areas:
    if area.type == 'VIEW_3D':
        area.spaces[0].region_3d.view_perspective = 'CAMERA'
        break

# レンダリングプレビュー出力
output_dir = os.path.join(os.path.dirname(bpy.data.filepath), '')
if not output_dir:
    output_dir = '/Users/naoki/dev/lab/love/ver46/results/'

bpy.context.scene.render.filepath = os.path.join(output_dir, 'blast_preview.png')
bpy.context.scene.render.resolution_x = 960
bpy.context.scene.render.resolution_y = 540
bpy.context.scene.frame_set(3)
bpy.ops.render.render(write_still=True)

print(f"\nプレビュー保存: {bpy.context.scene.render.filepath}")
print("\n手動での操作:")
print("  1. ビューポート右上の丸いアイコン群の3番目(球体)をクリック → マテリアルプレビュー")
print("  2. テンキー0 → カメラビュー")
print("  3. スペースキー → アニメーション再生")
