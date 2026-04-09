#!/usr/bin/env python3
"""
爆風シミュレーション結果を ParaView 用 VTK ファイルにエクスポート

1D球対称データを3Dカーテシアングリッドに展開し、
各タイムステップの .vtr (VTK Rectilinear Grid) を出力。
+ .pvd (ParaView Data) タイムシリーズファイル。

Usage: python3 export_paraview.py
Output: results/vtk/blast_*.vtr + blast.pvd

ParaViewで開く: File → Open → blast.pvd
"""

import numpy as np
import os
import vtk
from vtk.util.numpy_support import numpy_to_vtk
from blast_simulation import BlastParams, EulerSolver, SedovTaylor


def run_simulation():
    """シミュレーション実行 (細かい時間ステップ)"""
    params = BlastParams(
        E0=4.184e9, rho0=1.225, p0=1.01325e5, gamma=1.4,
        r_max=150.0, Nr=3000, t_end=0.06, CFL=0.45,
    )

    print("シミュレーション実行中...")
    solver = EulerSolver(params)
    # 20ステップ (ParaViewアニメーション用)
    t_snap = np.linspace(0.5e-3, params.t_end, 20)
    snapshots = solver.run(t_snap)
    return snapshots, params


def create_3d_grid(snapshot, params, N=100):
    """
    1D球対称データを3Dカーテシアングリッドに展開

    N: 各軸の格子点数 (N x N x N/2 — 上半分のみ)
    """
    r_1d = snapshot['r']
    rho_1d = snapshot['rho']
    u_1d = snapshot['u']
    p_1d = snapshot['p']

    # 3Dグリッド範囲
    L = params.r_max * 0.5  # 表示範囲
    Nz = N // 2  # z >= 0 のみ (地面反射)

    x = np.linspace(-L, L, N)
    y = np.linspace(-L, L, N)
    z = np.linspace(0, L, Nz)

    # 各格子点の原点からの距離
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)

    # 1Dデータを3Dに補間
    rho_3d = np.interp(R.ravel(), r_1d, rho_1d).reshape(R.shape)
    p_3d = np.interp(R.ravel(), r_1d, p_1d).reshape(R.shape)
    u_3d = np.interp(R.ravel(), r_1d, u_1d).reshape(R.shape)

    # 過圧
    overpressure = p_3d - params.p0

    # 速度ベクトル (半径方向)
    R_safe = np.maximum(R, 1e-10)
    ux = u_3d * X / R_safe
    uy = u_3d * Y / R_safe
    uz = u_3d * Z / R_safe

    return {
        'x': x, 'y': y, 'z': z,
        'density': rho_3d,
        'pressure': p_3d,
        'overpressure': overpressure,
        'velocity_magnitude': u_3d,
        'velocity_x': ux,
        'velocity_y': uy,
        'velocity_z': uz,
        'log_pressure': np.log10(np.maximum(p_3d, 1.0)),
    }


def write_vtr(grid_data, filepath):
    """VTK Rectilinear Grid (.vtr) ファイルを出力"""
    x = grid_data['x']
    y = grid_data['y']
    z = grid_data['z']
    Nx, Ny, Nz = len(x), len(y), len(z)

    # VTK Rectilinear Grid
    rgrid = vtk.vtkRectilinearGrid()
    rgrid.SetDimensions(Nx, Ny, Nz)

    xcoords = vtk.vtkFloatArray()
    for v in x:
        xcoords.InsertNextValue(v)
    ycoords = vtk.vtkFloatArray()
    for v in y:
        ycoords.InsertNextValue(v)
    zcoords = vtk.vtkFloatArray()
    for v in z:
        zcoords.InsertNextValue(v)

    rgrid.SetXCoordinates(xcoords)
    rgrid.SetYCoordinates(ycoords)
    rgrid.SetZCoordinates(zcoords)

    # スカラーデータを追加
    scalar_fields = [
        ('Density', 'density'),
        ('Pressure', 'pressure'),
        ('Overpressure', 'overpressure'),
        ('VelocityMagnitude', 'velocity_magnitude'),
        ('LogPressure', 'log_pressure'),
    ]

    for vtk_name, key in scalar_fields:
        arr = numpy_to_vtk(grid_data[key].ravel(order='F'), deep=True)
        arr.SetName(vtk_name)
        rgrid.GetPointData().AddArray(arr)

    # 速度ベクトル
    vel = np.zeros((Nx * Ny * Nz, 3))
    vel[:, 0] = grid_data['velocity_x'].ravel(order='F')
    vel[:, 1] = grid_data['velocity_y'].ravel(order='F')
    vel[:, 2] = grid_data['velocity_z'].ravel(order='F')
    vel_vtk = numpy_to_vtk(vel, deep=True)
    vel_vtk.SetName('Velocity')
    vel_vtk.SetNumberOfComponents(3)
    rgrid.GetPointData().AddArray(vel_vtk)

    # 最初のスカラーをアクティブに
    rgrid.GetPointData().SetActiveScalars('Pressure')

    # ファイル出力
    writer = vtk.vtkXMLRectilinearGridWriter()
    writer.SetFileName(filepath)
    writer.SetInputData(rgrid)
    writer.SetDataModeToBinary()
    writer.Write()


def write_pvd(times, vtr_files, filepath):
    """ParaView Data (.pvd) タイムシリーズファイルを出力"""
    lines = ['<?xml version="1.0"?>',
             '<VTKFile type="Collection" version="0.1">',
             '  <Collection>']
    for t, f in zip(times, vtr_files):
        fname = os.path.basename(f)
        lines.append(f'    <DataSet timestep="{t:.6e}" file="vtk/{fname}"/>')
    lines.append('  </Collection>')
    lines.append('</VTKFile>')

    with open(filepath, 'w') as f:
        f.write('\n'.join(lines))


def write_paraview_state(pvd_path, state_path, params):
    """ParaView Python スクリプト (自動可視化)"""
    script = f'''#!/usr/bin/env pvpython
"""
ParaView 可視化スクリプト
Usage: pvpython paraview_visualize.py
  または ParaView の Python Shell で実行
"""
from paraview.simple import *

# データ読み込み
reader = OpenDataFile("{pvd_path}")
UpdatePipeline()

# カラーマップ設定 (圧力)
display = Show(reader)
ColorBy(display, ('POINTS', 'Pressure'))
pressureLUT = GetColorTransferFunction('Pressure')
pressureLUT.ApplyPreset('Cool to Warm', True)
pressureLUT.RescaleTransferFunction({params.p0}, {params.p0 * 50})

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
'''
    with open(state_path, 'w') as f:
        f.write(script)
    os.chmod(state_path, 0o755)


def main():
    output_dir = 'results'
    vtk_dir = os.path.join(output_dir, 'vtk')
    os.makedirs(vtk_dir, exist_ok=True)

    # シミュレーション
    snapshots, params = run_simulation()

    # 3D VTK エクスポート
    N = 80  # 80x80x40 格子 (~256K点)
    print(f"\n3D VTK エクスポート ({N}x{N}x{N//2} 格子)...")

    times = []
    vtr_files = []

    for i, snap in enumerate(snapshots):
        t = snap['t']
        times.append(t)

        print(f"  [{i+1}/{len(snapshots)}] t={t*1e3:.2f} ms → 3D変換中...", end='')
        grid = create_3d_grid(snap, params, N=N)

        vtr_path = os.path.join(vtk_dir, f'blast_{i:04d}.vtr')
        write_vtr(grid, vtr_path)
        vtr_files.append(vtr_path)
        print(f" 保存: {os.path.basename(vtr_path)}")

    # PVD (タイムシリーズ)
    pvd_path = os.path.join(output_dir, 'blast.pvd')
    write_pvd(times, vtr_files, pvd_path)
    print(f"\nPVD: {pvd_path}")

    # ParaView 自動化スクリプト
    pv_script = os.path.join(output_dir, 'paraview_visualize.py')
    write_paraview_state(
        os.path.abspath(pvd_path), pv_script, params)
    print(f"ParaViewスクリプト: {pv_script}")

    # ファイルサイズ
    total_mb = sum(os.path.getsize(f) for f in vtr_files) / 1e6
    print(f"\n合計データサイズ: {total_mb:.1f} MB")

    print(f"\n{'='*60}")
    print("ParaView で開く方法:")
    print(f"{'='*60}")
    print(f"  方法1: ParaView GUI")
    print(f"    File → Open → {os.path.abspath(pvd_path)}")
    print(f"    Apply → 再生ボタン")
    print(f"")
    print(f"  方法2: コマンドライン")
    print(f"    paraview {pvd_path}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
