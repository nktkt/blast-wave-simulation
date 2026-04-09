#!/usr/bin/env python3
"""
Blender 爆風可視化スクリプト

同心球殻を使った爆風シミュレーションの3D可視化。
各シェルの色と透明度を圧力データに基づいてアニメーション。

========================================
使い方:
========================================

方法1: コマンドライン (バックグラウンドレンダリング)
  blender --background --python blender_blast.py

方法2: コマンドライン (GUIで確認)
  blender --python blender_blast.py

方法3: Blender GUI内
  1. Blender を開く
  2. Scripting ワークスペースに切り替え
  3. このスクリプトを開いて実行 (▶ボタン)

========================================
事前準備:
========================================
  python3 export_blender_data.py
  → results/blast_blender.npz が生成される
"""

import bpy
import bmesh
import numpy as np
import os
import math

# =============================================================================
# 設定
# =============================================================================

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'results', 'blast_blender.npz')
SHELL_SUBDIVISIONS = 3      # Icosphere の分割数 (3=320面, 4=1280面)
SHOCK_SUBDIVISIONS = 4
SHELL_BASE_ALPHA = 0.04     # 各シェルの基本透明度
EMISSION_STRENGTH = 5.0     # 発光強度
RENDER_ENGINE = 'BLENDER_EEVEE'  # 'BLENDER_EEVEE' or 'CYCLES'
RENDER_SAMPLES = 64
RESOLUTION_X = 1920
RESOLUTION_Y = 1080


# =============================================================================
# データ読み込み
# =============================================================================

def load_data(path):
    """blast_blender.npz を読み込み"""
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} が見つかりません。\n"
            f"先に python3 export_blender_data.py を実行してください。"
        )
    data = dict(np.load(path))
    print(f"データ読み込み: {path}")
    print(f"  シェル: {len(data['r_shells'])}, "
          f"フレーム: {int(data['n_frames'])}")
    return data


# =============================================================================
# シーンセットアップ
# =============================================================================

def clear_scene():
    """既存オブジェクトをすべて削除"""
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)

    # 既存マテリアルも削除
    for mat in bpy.data.materials:
        bpy.data.materials.remove(mat)
    for mesh in bpy.data.meshes:
        bpy.data.meshes.remove(mesh)


def setup_scene(data):
    """シーン全体の設定"""
    scene = bpy.context.scene
    n_frames = int(data['n_frames'])
    fps = int(data['fps'])

    # フレーム設定
    scene.frame_start = 1
    scene.frame_end = n_frames
    scene.frame_current = 1
    scene.render.fps = fps

    # レンダリング設定
    scene.render.engine = RENDER_ENGINE
    scene.render.resolution_x = RESOLUTION_X
    scene.render.resolution_y = RESOLUTION_Y

    if RENDER_ENGINE == 'CYCLES':
        scene.cycles.samples = RENDER_SAMPLES
        scene.cycles.use_denoising = True
    else:
        # EEVEE 設定 (バージョン互換)
        eevee = scene.eevee
        if hasattr(eevee, 'use_bloom'):
            eevee.use_bloom = True
            eevee.bloom_threshold = 0.8
            eevee.bloom_intensity = 0.5

    # ワールド背景 (暗い)
    world = bpy.data.worlds.get('World')
    if world is None:
        world = bpy.data.worlds.new('World')
    scene.world = world
    world.use_nodes = True
    bg_node = world.node_tree.nodes.get('Background')
    if bg_node:
        bg_node.inputs['Color'].default_value = (0.01, 0.01, 0.02, 1)
        bg_node.inputs['Strength'].default_value = 0.1

    # 地面 (参照用)
    bpy.ops.mesh.primitive_plane_add(size=300, location=(0, 0, -0.1))
    ground = bpy.context.active_object
    ground.name = 'Ground'
    mat_ground = bpy.data.materials.new('GroundMat')
    mat_ground.use_nodes = True
    bsdf = mat_ground.node_tree.nodes.get('Principled BSDF')
    if bsdf:
        bsdf.inputs['Base Color'].default_value = (0.15, 0.15, 0.15, 1)
        bsdf.inputs['Roughness'].default_value = 0.9
    ground.data.materials.append(mat_ground)


def setup_camera(data):
    """カメラ配置"""
    r_max = float(np.max(data['shock_radii'])) * 2.5
    cam_dist = max(r_max, 60)

    bpy.ops.object.camera_add(
        location=(cam_dist * 0.8, -cam_dist * 0.9, cam_dist * 0.35)
    )
    cam = bpy.context.active_object
    cam.name = 'BlastCamera'

    # 原点を向く
    constraint = cam.constraints.new(type='TRACK_TO')
    constraint.target = _get_or_create_empty('Origin', (0, 0, 5))
    constraint.track_axis = 'TRACK_NEGATIVE_Z'
    constraint.up_axis = 'UP_Y'

    bpy.context.scene.camera = cam

    # カメラ設定
    cam.data.lens = 35
    cam.data.clip_end = 1000


def setup_lights(data):
    """照明"""
    r = float(np.max(data['shock_radii'])) * 2

    # キーライト (上方)
    bpy.ops.object.light_add(type='AREA', location=(r, -r, r * 1.5))
    light = bpy.context.active_object
    light.name = 'KeyLight'
    light.data.energy = 5000
    light.data.size = r

    # フィルライト
    bpy.ops.object.light_add(type='AREA', location=(-r * 0.8, r * 0.5, r * 0.3))
    light2 = bpy.context.active_object
    light2.name = 'FillLight'
    light2.data.energy = 1000
    light2.data.size = r * 0.5


def _get_or_create_empty(name, location=(0, 0, 0)):
    """空のオブジェクトを取得または作成"""
    obj = bpy.data.objects.get(name)
    if obj is None:
        bpy.ops.object.empty_add(location=location)
        obj = bpy.context.active_object
        obj.name = name
    return obj


# =============================================================================
# シェル (同心球殻) の作成
# =============================================================================

def create_shell_material(name):
    """
    圧力マッピング用マテリアルを作成

    構造:
      Value (keyframed) → ColorRamp → Emission + Transparent → Mix → Output
    """
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    if hasattr(mat, 'blend_method'):
        mat.blend_method = 'BLEND'
    if hasattr(mat, 'shadow_method'):
        mat.shadow_method = 'NONE'
    mat.use_backface_culling = True

    tree = mat.node_tree
    tree.nodes.clear()

    # ノード作成
    output = tree.nodes.new('ShaderNodeOutputMaterial')
    output.location = (600, 0)

    mix = tree.nodes.new('ShaderNodeMixShader')
    mix.location = (400, 0)

    transparent = tree.nodes.new('ShaderNodeBsdfTransparent')
    transparent.location = (200, 100)

    emission = tree.nodes.new('ShaderNodeEmission')
    emission.location = (200, -100)
    emission.inputs['Strength'].default_value = EMISSION_STRENGTH

    ramp = tree.nodes.new('ShaderNodeValToRGB')
    ramp.location = (0, -100)

    # カラーランプ設定 (青→橙→白)
    cr = ramp.color_ramp
    cr.elements[0].position = 0.0
    cr.elements[0].color = (0.0, 0.0, 0.1, 1)   # 暗い青 (低圧)
    cr.elements[1].position = 1.0
    cr.elements[1].color = (1.0, 1.0, 0.9, 1)    # 白 (高圧)

    e1 = cr.elements.new(0.2)
    e1.color = (0.0, 0.05, 0.4, 1)  # 青
    e2 = cr.elements.new(0.45)
    e2.color = (0.8, 0.2, 0.0, 1)   # 橙
    e3 = cr.elements.new(0.7)
    e3.color = (1.0, 0.7, 0.1, 1)   # 黄

    # Value ノード (アニメーション制御用)
    value_node = tree.nodes.new('ShaderNodeValue')
    value_node.location = (-200, -100)
    value_node.outputs[0].default_value = 0.0
    value_node.name = 'PressureValue'
    value_node.label = 'Pressure'

    # Alpha 制御用 Math ノード
    math_alpha = tree.nodes.new('ShaderNodeMath')
    math_alpha.location = (-200, 100)
    math_alpha.operation = 'MULTIPLY'
    math_alpha.inputs[1].default_value = SHELL_BASE_ALPHA * 10

    # リンク
    tree.links.new(value_node.outputs[0], ramp.inputs['Fac'])
    tree.links.new(value_node.outputs[0], math_alpha.inputs[0])
    tree.links.new(ramp.outputs['Color'], emission.inputs['Color'])
    tree.links.new(math_alpha.outputs[0], mix.inputs['Fac'])
    tree.links.new(transparent.outputs[0], mix.inputs[1])
    tree.links.new(emission.outputs[0], mix.inputs[2])
    tree.links.new(mix.outputs[0], output.inputs['Surface'])

    return mat


def create_shells(data):
    """同心球殻メッシュを作成"""
    r_shells = data['r_shells']
    n_shells = len(r_shells)
    print(f"シェル作成中: {n_shells} 個...")

    # 親オブジェクト
    parent = _get_or_create_empty('BlastShells', (0, 0, 0))

    shells = []
    for i, r in enumerate(r_shells):
        if r < 0.3:
            continue

        bpy.ops.mesh.primitive_ico_sphere_add(
            subdivisions=SHELL_SUBDIVISIONS,
            radius=r,
            location=(0, 0, 0),
        )
        obj = bpy.context.active_object
        obj.name = f'Shell_{i:03d}'
        obj.parent = parent

        # スムーズシェーディング
        bpy.ops.object.shade_smooth()

        # マテリアル
        mat = create_shell_material(f'ShellMat_{i:03d}')
        obj.data.materials.append(mat)

        shells.append(obj)

        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{n_shells} シェル作成完了")

    print(f"  全 {len(shells)} シェル作成完了")
    return shells


def create_shock_sphere(data):
    """衝撃波面の発光球を作成"""
    bpy.ops.mesh.primitive_ico_sphere_add(
        subdivisions=SHOCK_SUBDIVISIONS,
        radius=1.0,
        location=(0, 0, 0),
    )
    shock = bpy.context.active_object
    shock.name = 'ShockFront'
    bpy.ops.object.shade_smooth()

    # 発光マテリアル
    mat = bpy.data.materials.new('ShockMat')
    mat.use_nodes = True
    if hasattr(mat, 'blend_method'):
        mat.blend_method = 'BLEND'
    if hasattr(mat, 'shadow_method'):
        mat.shadow_method = 'NONE'

    tree = mat.node_tree
    tree.nodes.clear()

    output = tree.nodes.new('ShaderNodeOutputMaterial')
    output.location = (400, 0)

    mix = tree.nodes.new('ShaderNodeMixShader')
    mix.location = (200, 0)
    mix.inputs['Fac'].default_value = 0.15  # ほぼ透明

    transparent = tree.nodes.new('ShaderNodeBsdfTransparent')
    transparent.location = (0, 100)

    emission = tree.nodes.new('ShaderNodeEmission')
    emission.location = (0, -100)
    emission.inputs['Color'].default_value = (0.3, 0.7, 1.0, 1)  # 水色
    emission.inputs['Strength'].default_value = 3.0

    tree.links.new(transparent.outputs[0], mix.inputs[1])
    tree.links.new(emission.outputs[0], mix.inputs[2])
    tree.links.new(mix.outputs[0], output.inputs['Surface'])

    shock.data.materials.append(mat)
    return shock


# =============================================================================
# アニメーション
# =============================================================================

def animate_shells(shells, data):
    """各シェルの圧力値をフレームごとにキーフレーム設定"""
    p_norm = data['p_norm']
    n_frames = int(data['n_frames'])
    r_shells = data['r_shells']
    n_shells = len(r_shells)

    print(f"アニメーション設定中: {n_frames} フレーム x {len(shells)} シェル...")

    shell_idx_map = {}
    for obj in shells:
        idx = int(obj.name.split('_')[1])
        shell_idx_map[idx] = obj

    for frame in range(n_frames):
        frame_num = frame + 1  # Blender のフレームは1始まり

        for idx, obj in shell_idx_map.items():
            if idx >= p_norm.shape[1]:
                continue

            mat = obj.data.materials[0]
            value_node = mat.node_tree.nodes.get('PressureValue')
            if value_node is None:
                continue

            val = float(p_norm[frame, idx])
            value_node.outputs[0].default_value = val
            value_node.outputs[0].keyframe_insert(
                data_path='default_value', frame=frame_num
            )

        if (frame + 1) % 10 == 0:
            print(f"  フレーム {frame+1}/{n_frames} 完了")

    # キーフレーム補間を LINEAR に設定 (Blender 5.x 互換)
    try:
        for obj in shells:
            mat = obj.data.materials[0]
            anim = mat.node_tree.animation_data
            if anim and anim.action:
                action = anim.action
                # Blender 5.x: layered action system
                if hasattr(action, 'layers'):
                    for layer in action.layers:
                        for strip in layer.strips:
                            for cb in strip.channelbags:
                                for fc in cb.fcurves:
                                    for kp in fc.keyframe_points:
                                        kp.interpolation = 'LINEAR'
                # Blender 4.x 以前
                elif hasattr(action, 'fcurves'):
                    for fc in action.fcurves:
                        for kp in fc.keyframe_points:
                            kp.interpolation = 'LINEAR'
    except Exception as e:
        print(f"  補間設定スキップ: {e}")

    print("  アニメーション設定完了")


def animate_shock_sphere(shock_obj, data):
    """衝撃波面球のスケールアニメーション"""
    shock_radii = data['shock_radii']
    n_frames = int(data['n_frames'])

    print("衝撃波面アニメーション設定中...")

    for frame in range(n_frames):
        frame_num = frame + 1
        r = float(shock_radii[frame])
        shock_obj.scale = (r, r, r)
        shock_obj.keyframe_insert(data_path='scale', frame=frame_num)

    # LINEAR 補間 (Blender 5.x 互換)
    try:
        anim = shock_obj.animation_data
        if anim and anim.action:
            action = anim.action
            if hasattr(action, 'layers'):
                for layer in action.layers:
                    for strip in layer.strips:
                        for cb in strip.channelbags:
                            for fc in cb.fcurves:
                                for kp in fc.keyframe_points:
                                    kp.interpolation = 'LINEAR'
            elif hasattr(action, 'fcurves'):
                for fc in action.fcurves:
                    for kp in fc.keyframe_points:
                        kp.interpolation = 'LINEAR'
    except Exception as e:
        print(f"  補間設定スキップ: {e}")

    print("  衝撃波面アニメーション完了")


# =============================================================================
# テキスト情報
# =============================================================================

def add_info_text(data):
    """画面にパラメータ情報テキストを追加"""
    E0 = float(data['E0'])
    W_kg = E0 / 4.184e6

    bpy.ops.object.text_add(location=(-30, -50, 0.1))
    text_obj = bpy.context.active_object
    text_obj.name = 'InfoText'
    text_obj.data.body = (
        f"Blast Wave Simulation\n"
        f"E0 = {E0:.2e} J ({W_kg:.0f} kg TNT)\n"
        f"Sedov-Taylor Solution"
    )
    text_obj.data.size = 2.0
    text_obj.rotation_euler = (math.radians(90), 0, 0)

    mat = bpy.data.materials.new('TextMat')
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get('Principled BSDF')
    if bsdf:
        bsdf.inputs['Base Color'].default_value = (0.8, 0.8, 0.8, 1)
        bsdf.inputs['Emission Color'].default_value = (0.8, 0.8, 0.8, 1)
        bsdf.inputs['Emission Strength'].default_value = 1.0
    text_obj.data.materials.append(mat)


# =============================================================================
# メイン
# =============================================================================

def main():
    print("=" * 60)
    print("Blender 爆風可視化")
    print("=" * 60)

    # データ読み込み
    data = load_data(DATA_PATH)

    # シーン構築
    print("\nシーン構築中...")
    clear_scene()
    setup_scene(data)
    setup_camera(data)
    setup_lights(data)

    # オブジェクト作成
    shells = create_shells(data)
    shock = create_shock_sphere(data)
    add_info_text(data)

    # アニメーション
    print("\nアニメーション設定中...")
    animate_shells(shells, data)
    animate_shock_sphere(shock, data)

    # フレーム1に戻る
    bpy.context.scene.frame_set(1)

    # .blend ファイル保存
    blend_path = os.path.join(os.path.dirname(DATA_PATH), 'blast_wave.blend')
    bpy.ops.wm.save_as_mainfile(filepath=blend_path)
    print(f"\n保存: {blend_path}")

    print("\n" + "=" * 60)
    print("完了！")
    print("=" * 60)
    print(f"\n操作方法:")
    print(f"  - スペースキー: アニメーション再生")
    print(f"  - テンキー0: カメラビュー")
    print(f"  - F12: レンダリング")
    print(f"  - Ctrl+F12: アニメーションレンダリング")


if __name__ == '__main__':
    main()
