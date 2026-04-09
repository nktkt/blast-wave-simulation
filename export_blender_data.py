#!/usr/bin/env python3
"""
blast_simulation の結果を Blender 用に .npz エクスポート

Usage: python3 export_blender_data.py
Output: results/blast_blender.npz
"""

import numpy as np
import os
from blast_simulation import BlastParams, EulerSolver, SedovTaylor

N_SHELLS = 120      # Blenderに渡す半径方向のシェル数
N_FRAMES = 60       # アニメーションフレーム数
FPS = 24


def export():
    params = BlastParams(
        E0=4.184e9, rho0=1.225, p0=1.01325e5, gamma=1.4,
        r_max=150.0, Nr=3000, t_end=0.08, CFL=0.45,
    )

    # 細かい時間ステップでシミュレーション (スムーズなアニメーション用)
    print(f"シミュレーション実行 ({N_FRAMES} フレーム)...")
    solver = EulerSolver(params)
    t_snap = np.linspace(0.5e-3, params.t_end, N_FRAMES)
    snapshots = solver.run(t_snap)

    # 半径方向をダウンサンプル (衝撃波付近を密に)
    r_full = snapshots[0]['r']
    r_max_vis = params.r_max * 0.5
    r_shells = np.concatenate([
        np.linspace(0.5, 40, N_SHELLS * 2 // 3),
        np.linspace(40, r_max_vis, N_SHELLS // 3 + 1)[1:],
    ])[:N_SHELLS]

    times = np.array([s['t'] for s in snapshots])

    # 各フレームの物理量を補間
    p_data = np.zeros((N_FRAMES, N_SHELLS))
    rho_data = np.zeros((N_FRAMES, N_SHELLS))
    u_data = np.zeros((N_FRAMES, N_SHELLS))
    shock_radii = np.zeros(N_FRAMES)

    sedov = SedovTaylor(params)

    for i, s in enumerate(snapshots):
        p_data[i] = np.interp(r_shells, s['r'], s['p'])
        rho_data[i] = np.interp(r_shells, s['r'], s['rho'])
        u_data[i] = np.interp(r_shells, s['r'], s['u'])

        # 衝撃波位置
        above = s['p'] > params.p0 * 1.5
        if np.any(above):
            shock_radii[i] = s['r'][np.where(above)[0][-1]]
        else:
            shock_radii[i] = sedov.shock_radius(s['t'])

    # 正規化 (log スケール, 0~1)
    p_log = np.log10(np.maximum(p_data / params.p0, 1.0))  # log10(p/p0), >=0
    p_log_max = max(np.max(p_log), 1.0)
    p_norm = np.clip(p_log / p_log_max, 0, 1)

    rho_norm = np.clip((rho_data - params.rho0) / (6 * params.rho0), 0, 1)

    output_path = os.path.join('results', 'blast_blender.npz')
    np.savez(
        output_path,
        r_shells=r_shells,
        times=times,
        p_norm=p_norm,
        rho_norm=rho_norm,
        u_data=u_data,
        shock_radii=shock_radii,
        p_log_max=p_log_max,
        E0=params.E0,
        p0=params.p0,
        rho0=params.rho0,
        gamma=params.gamma,
        n_frames=N_FRAMES,
        fps=FPS,
    )

    print(f"\n保存: {output_path}")
    print(f"  シェル数: {N_SHELLS}, フレーム数: {N_FRAMES}, FPS: {FPS}")
    print(f"  r_shells: {r_shells[0]:.1f} ~ {r_shells[-1]:.1f} m")
    print(f"  衝撃波半径: {shock_radii[0]:.1f} ~ {np.max(shock_radii):.1f} m")
    print(f"  p_log_max: {p_log_max:.2f}")
    print(f"\n次のステップ:")
    print(f"  blender --python blender_blast.py")
    print(f"  または Blender を開いて Scripting タブで blender_blast.py を実行")


if __name__ == '__main__':
    export()
