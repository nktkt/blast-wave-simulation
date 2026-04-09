#!/usr/bin/env python3
"""
爆風可視化 v2 — 鮮やかな色付きの爆風アニメーション

使い方:
  /Applications/Blender.app/Contents/MacOS/Blender --python blender_blast_v2.py
"""

import bpy
import numpy as np
import os
import math

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'results', 'blast_blender.npz')


def clear_scene():
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()
    for block in [bpy.data.materials, bpy.data.meshes]:
        for item in block:
            block.remove(item)


def load_data():
    data = dict(np.load(DATA_PATH))
    print(f"データ: {int(data['n_frames'])} フレーム, "
          f"衝撃波 {data['shock_radii'][0]:.1f}~{np.max(data['shock_radii']):.1f} m")
    return data


def make_emission_material(name, color_rgb, strength=5.0):
    """発光マテリアル (レンダリングで色が出る)"""
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = (*color_rgb, 1.0)  # ソリッドモード用
    tree = mat.node_tree
    nodes = tree.nodes
    links = tree.links
    nodes.clear()

    output = nodes.new('ShaderNodeOutputMaterial')
    output.location = (300, 0)

    emission = nodes.new('ShaderNodeEmission')
    emission.location = (0, 0)
    emission.inputs['Color'].default_value = (*color_rgb, 1.0)
    emission.inputs['Strength'].default_value = strength

    links.new(emission.outputs[0], output.inputs['Surface'])
    return mat


def make_transparent_emission_material(name, color_rgb, strength=3.0, alpha=0.6):
    """半透明発光マテリアル"""
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = (*color_rgb, alpha)
    tree = mat.node_tree
    nodes = tree.nodes
    links = tree.links
    nodes.clear()

    output = nodes.new('ShaderNodeOutputMaterial')
    output.location = (500, 0)

    mix = nodes.new('ShaderNodeMixShader')
    mix.location = (300, 0)
    mix.inputs['Fac'].default_value = alpha

    transparent = nodes.new('ShaderNodeBsdfTransparent')
    transparent.location = (100, 100)

    emission = nodes.new('ShaderNodeEmission')
    emission.location = (100, -100)
    emission.inputs['Color'].default_value = (*color_rgb, 1.0)
    emission.inputs['Strength'].default_value = strength

    links.new(transparent.outputs[0], mix.inputs[1])
    links.new(emission.outputs[0], mix.inputs[2])
    links.new(mix.outputs[0], output.inputs['Surface'])

    if hasattr(mat, 'blend_method'):
        mat.blend_method = 'BLEND'
    return mat


def setup_world():
    scene = bpy.context.scene
    scene.render.engine = 'BLENDER_EEVEE'
    scene.render.resolution_x = 1920
    scene.render.resolution_y = 1080
    scene.render.fps = 24

    world = bpy.data.worlds.get('World') or bpy.data.worlds.new('World')
    scene.world = world
    tree = world.node_tree
    bg = tree.nodes.get('Background')
    if bg:
        bg.inputs['Color'].default_value = (0.02, 0.02, 0.04, 1)
        bg.inputs['Strength'].default_value = 1.0


def create_ground():
    bpy.ops.mesh.primitive_plane_add(size=400, location=(0, 0, -0.05))
    ground = bpy.context.active_object
    ground.name = 'Ground'
    mat = bpy.data.materials.new('GroundMat')
    mat.diffuse_color = (0.08, 0.08, 0.08, 1)
    tree = mat.node_tree
    nodes = tree.nodes
    bsdf = nodes.get('Principled BSDF')
    if bsdf:
        bsdf.inputs['Base Color'].default_value = (0.08, 0.08, 0.08, 1)
        bsdf.inputs['Roughness'].default_value = 0.95
    ground.data.materials.append(mat)
    ground.color = (0.08, 0.08, 0.08, 1)


def create_shock_wave(data):
    """膨張する衝撃波球 + 高圧コア"""
    n_frames = int(data['n_frames'])
    shock_radii = data['shock_radii']

    # --- 衝撃波面 (オレンジ発光) ---
    bpy.ops.mesh.primitive_uv_sphere_add(
        segments=64, ring_count=32, radius=1.0, location=(0, 0, 0))
    shock = bpy.context.active_object
    shock.name = 'ShockWave'
    bpy.ops.object.shade_smooth()
    shock.color = (1.0, 0.4, 0.05, 1.0)

    mat_shock = make_transparent_emission_material(
        'ShockMat', (1.0, 0.4, 0.05), strength=8.0, alpha=0.4)
    shock.data.materials.append(mat_shock)

    for f in range(n_frames):
        r = max(float(shock_radii[f]), 0.1)
        shock.scale = (r, r, r)
        shock.keyframe_insert(data_path='scale', frame=f + 1)

    # --- 高圧コア (白~黄 発光) ---
    bpy.ops.mesh.primitive_uv_sphere_add(
        segments=48, ring_count=24, radius=1.0, location=(0, 0, 0))
    core = bpy.context.active_object
    core.name = 'Core'
    bpy.ops.object.shade_smooth()
    core.color = (1.0, 0.95, 0.5, 1.0)

    mat_core = make_emission_material('CoreMat', (1.0, 0.95, 0.5), strength=15.0)
    core.data.materials.append(mat_core)

    for f in range(n_frames):
        r = max(float(shock_radii[f]) * 0.35, 0.05)
        core.scale = (r, r, r)
        core.keyframe_insert(data_path='scale', frame=f + 1)

    # --- 中間層 (赤~橙) ---
    bpy.ops.mesh.primitive_uv_sphere_add(
        segments=48, ring_count=24, radius=1.0, location=(0, 0, 0))
    mid = bpy.context.active_object
    mid.name = 'MidLayer'
    bpy.ops.object.shade_smooth()
    mid.color = (1.0, 0.15, 0.0, 1.0)

    mat_mid = make_transparent_emission_material(
        'MidMat', (1.0, 0.15, 0.0), strength=5.0, alpha=0.5)
    mid.data.materials.append(mat_mid)

    for f in range(n_frames):
        r = max(float(shock_radii[f]) * 0.7, 0.08)
        mid.scale = (r, r, r)
        mid.keyframe_insert(data_path='scale', frame=f + 1)


def create_ground_rings(data):
    """地面上の圧力波リング"""
    shock_radii = data['shock_radii']
    n_frames = int(data['n_frames'])

    ring_specs = [
        ('Ring_Shock', (1.0, 0.3, 0.0), 1.0, 6.0),
        ('Ring_Mid', (1.0, 0.6, 0.0), 1.4, 3.0),
        ('Ring_Low', (0.8, 0.8, 0.0), 1.8, 1.5),
    ]

    for name, color, scale_mult, emission_str in ring_specs:
        bpy.ops.mesh.primitive_torus_add(
            major_radius=1.0, minor_radius=0.15,
            major_segments=96, minor_segments=12,
            location=(0, 0, 0.1))
        ring = bpy.context.active_object
        ring.name = name
        ring.color = (*color, 1.0)

        mat = make_emission_material(f'{name}Mat', color, emission_str)
        ring.data.materials.append(mat)

        for f in range(n_frames):
            r = max(float(shock_radii[f]) * scale_mult, 0.1)
            ring.scale = (r, r, 0.5)
            ring.keyframe_insert(data_path='scale', frame=f + 1)


def create_distance_markers():
    """距離マーカー"""
    for dist in [10, 20, 30, 50]:
        bpy.ops.mesh.primitive_torus_add(
            major_radius=dist, minor_radius=0.05,
            major_segments=64, minor_segments=8,
            location=(0, 0, 0.02))
        marker = bpy.context.active_object
        marker.name = f'Dist_{dist}m'
        mat = make_emission_material(f'DistMat_{dist}', (0.3, 0.3, 0.4), 0.5)
        marker.data.materials.append(mat)
        marker.color = (0.3, 0.3, 0.4, 1)

        bpy.ops.object.text_add(location=(dist + 1, 0, 0.3))
        txt = bpy.context.active_object
        txt.name = f'Text_{dist}m'
        txt.data.body = f'{dist}m'
        txt.data.size = 1.5
        txt.rotation_euler = (math.radians(90), 0, 0)
        mat_t = make_emission_material(f'TextMat_{dist}', (0.6, 0.6, 0.7), 1.0)
        txt.data.materials.append(mat_t)
        txt.color = (0.6, 0.6, 0.7, 1)


def create_camera_and_lights(data):
    max_r = float(np.max(data['shock_radii']))
    cam_dist = max_r * 2.5

    bpy.ops.object.camera_add(
        location=(cam_dist * 0.5, -cam_dist * 0.7, cam_dist * 0.4))
    cam = bpy.context.active_object
    cam.name = 'Camera'
    cam.data.lens = 28
    cam.data.clip_end = 2000

    bpy.ops.object.empty_add(location=(0, 0, max_r * 0.25))
    target = bpy.context.active_object
    target.name = 'CamTarget'

    c = cam.constraints.new(type='TRACK_TO')
    c.target = target
    c.track_axis = 'TRACK_NEGATIVE_Z'
    c.up_axis = 'UP_Y'
    bpy.context.scene.camera = cam

    # 爆心のポイントライト (オレンジ)
    bpy.ops.object.light_add(type='POINT', location=(0, 0, 3))
    pl = bpy.context.active_object
    pl.name = 'ExplosionLight'
    pl.data.energy = 100000
    pl.data.color = (1.0, 0.5, 0.1)
    pl.data.shadow_soft_size = 10

    # 環境光 (弱い)
    bpy.ops.object.light_add(type='SUN', location=(30, -20, 50))
    sun = bpy.context.active_object
    sun.name = 'Sun'
    sun.data.energy = 0.5
    sun.data.color = (0.7, 0.8, 1.0)


def setup_viewport():
    """ビューポートをオブジェクトカラー表示に"""
    for area in bpy.context.screen.areas:
        if area.type == 'VIEW_3D':
            for space in area.spaces:
                if space.type == 'VIEW_3D':
                    space.shading.type = 'SOLID'
                    space.shading.color_type = 'OBJECT'
                    space.shading.light = 'FLAT'
                    space.clip_end = 2000
                    space.region_3d.view_perspective = 'CAMERA'


def render_frames(data):
    """プレビュー画像"""
    output_dir = os.path.dirname(DATA_PATH)
    frames = [1, 3, 8, 15, 30, 45, 60]
    frames = [f for f in frames if f <= int(data['n_frames'])]

    for f in frames:
        bpy.context.scene.frame_set(f)
        path = os.path.join(output_dir, f'blast_frame_{f:03d}.png')
        bpy.context.scene.render.filepath = path
        bpy.context.scene.render.resolution_x = 960
        bpy.context.scene.render.resolution_y = 540
        bpy.ops.render.render(write_still=True)
        print(f"  frame {f} → {path}")


def main():
    print("=" * 50)
    print("爆風可視化 v2")
    print("=" * 50)

    data = load_data()
    n = int(data['n_frames'])

    clear_scene()
    setup_world()
    bpy.context.scene.frame_start = 1
    bpy.context.scene.frame_end = n

    create_ground()
    create_shock_wave(data)
    create_ground_rings(data)
    create_distance_markers()
    create_camera_and_lights(data)

    try:
        setup_viewport()
    except Exception:
        pass

    bpy.context.scene.frame_set(5)

    blend_path = os.path.join(os.path.dirname(DATA_PATH), 'blast_wave_v2.blend')
    bpy.ops.wm.save_as_mainfile(filepath=blend_path)
    print(f"保存: {blend_path}")

    print("レンダリング中...")
    render_frames(data)
    print("完了！")


if __name__ == '__main__':
    main()
