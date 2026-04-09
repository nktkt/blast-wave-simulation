#!/usr/bin/env python3
"""
爆風シミュレーション (Blast Wave Simulation)
==============================================
球対称爆風の数値シミュレーションと Sedov-Taylor 解析解の比較

物理モデル:
  - 1D球対称 圧縮性オイラー方程式 (幾何学的ソース項)
  - 有限体積法 + HLL Riemann Solver
  - Sedov-Taylor 自己相似解 (解析解)
  - Rankine-Hugoniot 衝撃波関係式
  - Hopkinson-Cranz スケーリング則

Author: Research Simulation
Date: 2026-04-09
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os


# =============================================================================
# パラメータ
# =============================================================================

class BlastParams:
    def __init__(self, E0=4.184e9, rho0=1.225, p0=1.01325e5, gamma=1.4,
                 r_max=150.0, Nr=3000, t_end=0.08, CFL=0.45):
        self.E0 = E0
        self.rho0 = rho0
        self.p0 = p0
        self.gamma = gamma
        self.r_max = r_max
        self.Nr = Nr
        self.t_end = t_end
        self.CFL = CFL
        self.c0 = np.sqrt(gamma * p0 / rho0)


# =============================================================================
# Sedov-Taylor 解析解
# =============================================================================

class SedovTaylor:
    """Sedov-Taylor 自己相似解 (球対称)"""

    def __init__(self, params):
        self.p = params
        self.gamma = params.gamma
        # alpha: E = alpha * rho0 * R^5 / t^2
        g_tab = [1.2, 1.4, 5/3, 2.0, 3.0, 5.0, 7.0]
        a_tab = [1.033, 0.851, 0.493, 0.311, 0.153, 0.074, 0.050]
        self.alpha = float(np.interp(params.gamma, g_tab, a_tab))

    def shock_radius(self, t):
        return (self.p.E0 / (self.alpha * self.p.rho0)) ** 0.2 * t ** 0.4

    def shock_velocity(self, t):
        return 0.4 * self.shock_radius(t) / t

    def shock_mach(self, t):
        return self.shock_velocity(t) / self.p.c0

    def shock_overpressure(self, t):
        """衝撃波ピーク過圧 [Pa]"""
        g = self.gamma
        Ms = self.shock_mach(t)
        if Ms <= 1.0:
            return 0.0
        return self.p.p0 * (2 * g * Ms ** 2 - (g - 1)) / (g + 1) - self.p.p0


# =============================================================================
# 1D球対称 圧縮性オイラー方程式ソルバー
#
# 方程式系 (ソース項分離形式):
#   ∂U/∂t + ∂F/∂r = S
# ここで U=(ρ, ρu, E), F=(ρu, ρu²+p, (E+p)u)
# ソース: S = -(2/r)(ρu, ρu², (E+p)u)
#
# 手法: Strang splitting — フラックス (HLL) + ソース (半陰的)
# =============================================================================

class EulerSolver:
    def __init__(self, params):
        self.par = params
        self.g = params.gamma
        Nr = params.Nr

        # 等間隔格子
        self.dr = params.r_max / Nr
        self.r = (np.arange(Nr) + 0.5) * self.dr  # セル中心
        self.r_edge = np.arange(Nr + 1) * self.dr  # セル境界

        # 初期条件
        rho = np.full(Nr, params.rho0)
        u = np.zeros(Nr)
        p = np.full(Nr, params.p0)

        # 爆発エネルギーを小球に集中
        r_init = max(10.0 * self.dr, 0.5)
        vol_init = (4.0 / 3.0) * np.pi * r_init ** 3
        p_init = params.E0 * (params.gamma - 1) / vol_init
        mask = self.r <= r_init
        p[mask] = p_init

        print(f"初期条件: r_init={r_init:.3f} m, "
              f"p_init={p_init:.3e} Pa ({p_init/params.p0:.0f} atm), "
              f"dr={self.dr:.4f} m")

        # 保存変数 U = [rho, rho*u, E]  shape: (3, Nr)
        self.U = np.zeros((3, Nr))
        self.U[0] = rho
        self.U[1] = rho * u
        self.U[2] = p / (params.gamma - 1) + 0.5 * rho * u ** 2

    def _prim(self, U):
        """保存→プリミティブ"""
        rho = np.maximum(U[0], 1e-20)
        u = U[1] / rho
        p = (self.g - 1) * (U[2] - 0.5 * rho * u ** 2)
        p = np.maximum(p, 1e-10)
        return rho, u, p

    def _flux(self, rho, u, p):
        E = p / (self.g - 1) + 0.5 * rho * u ** 2
        return np.array([rho * u, rho * u ** 2 + p, (E + p) * u])

    def _cs(self, rho, p):
        return np.sqrt(self.g * p / rho)

    def _hll(self, rhoL, uL, pL, rhoR, uR, pR):
        """HLL Riemann solver"""
        cL = self._cs(rhoL, pL)
        cR = self._cs(rhoR, pR)
        SL = np.minimum(uL - cL, uR - cR)
        SR = np.maximum(uL + cL, uR + cR)

        FL = self._flux(rhoL, uL, pL)
        FR = self._flux(rhoR, uR, pR)
        UL = np.array([rhoL, rhoL * uL,
                        pL / (self.g - 1) + 0.5 * rhoL * uL ** 2])
        UR = np.array([rhoR, rhoR * uR,
                        pR / (self.g - 1) + 0.5 * rhoR * uR ** 2])

        denom = SR - SL + 1e-30
        F_hll = (SR * FL - SL * FR + SL * SR * (UR - UL)) / denom

        # 条件分岐
        mask_L = SL >= 0
        mask_R = SR <= 0
        F = np.where(mask_L, FL, np.where(mask_R, FR, F_hll))
        return F

    def _compute_dt(self):
        rho, u, p = self._prim(self.U)
        c = self._cs(rho, p)
        return self.par.CFL * self.dr / np.max(np.abs(u) + c)

    def _step(self, dt):
        """Strang splitting: S(dt/2) → F(dt) → S(dt/2)"""
        # --- ソース半ステップ ---
        self._apply_source(dt * 0.5)

        # --- フラックスステップ ---
        rho, u, p = self._prim(self.U)

        # ゴースト (反射 + 外部固定)
        rho_g = np.concatenate([[rho[0]], rho, [rho[-1]]])
        u_g = np.concatenate([[-u[0]], u, [u[-1]]])
        p_g = np.concatenate([[p[0]], p, [p[-1]]])

        # セル境界フラックス (i+1/2 for i=0..Nr)
        F_edge = self._hll(
            rho_g[:-1], u_g[:-1], p_g[:-1],
            rho_g[1:], u_g[1:], p_g[1:]
        )  # shape (3, Nr+1)

        # 更新 ∂U/∂t = -(F_{i+1/2} - F_{i-1/2}) / dr
        self.U -= dt / self.dr * (F_edge[:, 1:] - F_edge[:, :-1])

        # --- ソース半ステップ ---
        self._apply_source(dt * 0.5)

        # 物理制約
        self.U[0] = np.maximum(self.U[0], 1e-20)
        rho, u, p = self._prim(self.U)
        p = np.maximum(p, self.par.p0 * 1e-6)
        self.U[2] = p / (self.g - 1) + 0.5 * rho * u ** 2

    def _apply_source(self, dt):
        """
        幾何学的ソース項 S = -(2/r)(ρu, ρu², (E+p)u)
        半陰的スキームで適用 (安定性のため)
        """
        rho, u, p = self._prim(self.U)
        r = self.r
        E = self.U[2]

        # r=0 付近のゼロ割り防止
        inv_r = np.where(r > self.dr, 2.0 / r, 0.0)

        # ソース項 (半陰的: u_new = u_old / (1 + 2*|u|*dt/r) の形で減衰を安定化)
        factor = 1.0 / (1.0 + inv_r * np.abs(u) * dt)

        self.U[0] -= dt * inv_r * rho * u * factor
        self.U[1] -= dt * inv_r * rho * u ** 2 * factor
        self.U[2] -= dt * inv_r * (E + p) * u * factor

    def run(self, t_snapshots):
        """シミュレーション実行"""
        t_snapshots = np.sort(t_snapshots)
        snapshots = []
        t = 0.0
        nstep = 0
        snap_idx = 0

        while snap_idx < len(t_snapshots):
            dt = self._compute_dt()
            # 次のスナップショット時刻に合わせる
            if t + dt >= t_snapshots[snap_idx]:
                dt = t_snapshots[snap_idx] - t
                if dt < 1e-15:
                    dt = self._compute_dt()

            self._step(dt)
            t += dt
            nstep += 1

            if snap_idx < len(t_snapshots) and t >= t_snapshots[snap_idx] - 1e-12:
                rho, u, p = self._prim(self.U)
                snapshots.append({
                    't': t,
                    'r': self.r.copy(),
                    'rho': rho.copy(),
                    'u': u.copy(),
                    'p': p.copy(),
                })
                r_s, p_s = self._detect_shock(rho, u, p)
                print(f"  t={t*1e3:7.2f} ms | R_shock≈{r_s:6.1f} m | "
                      f"p_shock={p_s:.3e} Pa | "
                      f"max|u|={np.max(np.abs(u)):.0f} m/s | n={nstep}")
                snap_idx += 1

        print(f"計算完了: 総ステップ数 = {nstep}")
        return snapshots

    def _detect_shock(self, rho, u, p):
        """衝撃波位置を圧力勾配ピークから検出"""
        dp = np.abs(np.diff(p))
        # 原点付近の高圧領域を避けるため、圧力が大気圧の2倍以上→以下に
        # 遷移する点を探す
        p_threshold = self.par.p0 * 2
        above = p > p_threshold
        if not np.any(above):
            return 0.0, self.par.p0

        # 最も外側で閾値を超える点
        i_last = np.where(above)[0][-1]
        i_shock = min(i_last + 1, len(self.r) - 1)
        return self.r[i_shock], p[min(i_last, len(p) - 1)]


# =============================================================================
# Rankine-Hugoniot
# =============================================================================

class RankineHugoniot:
    def __init__(self, gamma=1.4):
        self.g = gamma

    def from_mach(self, Ms):
        g = self.g
        Ms = max(Ms, 1.0001)
        p_ratio = (2 * g * Ms ** 2 - (g - 1)) / (g + 1)
        rho_ratio = ((g + 1) * Ms ** 2) / ((g - 1) * Ms ** 2 + 2)
        T_ratio = p_ratio / rho_ratio
        return {'Ms': Ms, 'p_ratio': p_ratio, 'rho_ratio': rho_ratio,
                'T_ratio': T_ratio}


# =============================================================================
# 可視化
# =============================================================================

def plot_profiles(snapshots, params, output_dir):
    """密度・速度・圧力の空間分布"""
    sedov = SedovTaylor(params)
    n = len(snapshots)
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, n))

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle(
        f'Blast Wave — E₀ = {params.E0:.2e} J ({params.E0/4.184e6:.0f} kg TNT)',
        fontsize=14, fontweight='bold')

    for i, s in enumerate(snapshots):
        c = colors[i]
        lbl = f't={s["t"]*1e3:.1f} ms'
        R_s = sedov.shock_radius(s['t'])

        axes[0].plot(s['r'], s['rho'], '-', color=c, lw=1.2, label=lbl)
        axes[0].axvline(R_s, color=c, ls=':', alpha=0.3)
        axes[1].plot(s['r'], s['u'], '-', color=c, lw=1.2, label=lbl)
        axes[1].axvline(R_s, color=c, ls=':', alpha=0.3)
        axes[2].plot(s['r'], s['p'] / 1e3, '-', color=c, lw=1.2, label=lbl)
        axes[2].axvline(R_s, color=c, ls=':', alpha=0.3)

    labels = ['Density [kg/m³]', 'Velocity [m/s]', 'Pressure [kPa]']
    titles = ['ρ(r)', 'u(r)', 'p(r)']
    for ax, yl, tt in zip(axes, labels, titles):
        ax.set_xlabel('r [m]')
        ax.set_ylabel(yl)
        ax.set_title(tt)
        ax.legend(fontsize=7, loc='upper right')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, params.r_max * 0.6)

    plt.tight_layout()
    path = os.path.join(output_dir, 'blast_profiles.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  {path}")


def plot_shock_trajectory(snapshots, params, output_dir):
    """衝撃波軌跡の比較"""
    sedov = SedovTaylor(params)

    # 数値解の衝撃波位置
    num_t, num_R, num_P = [], [], []
    for s in snapshots:
        p_thresh = params.p0 * 2
        above = s['p'] > p_thresh
        if np.any(above):
            idx = np.where(above)[0][-1]
            num_R.append(s['r'][idx])
            num_P.append(s['p'][idx])
        else:
            num_R.append(0)
            num_P.append(params.p0)
        num_t.append(s['t'])

    num_t = np.array(num_t)
    num_R = np.array(num_R)
    num_P = np.array(num_P)

    # 解析解
    t_a = np.linspace(max(num_t[0] * 0.3, 1e-4), num_t[-1] * 1.1, 300)
    R_a = np.array([sedov.shock_radius(t) for t in t_a])
    V_a = np.array([sedov.shock_velocity(t) for t in t_a])
    P_a = np.array([sedov.shock_overpressure(t) + params.p0 for t in t_a])

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    fig.suptitle('Shock Trajectory — Numerical vs Sedov-Taylor',
                 fontsize=14, fontweight='bold')

    # R(t)
    axes[0].plot(t_a * 1e3, R_a, 'b-', lw=2, label='Sedov-Taylor')
    axes[0].plot(num_t * 1e3, num_R, 'ro', ms=6, label='Numerical')
    axes[0].set_xlabel('t [ms]')
    axes[0].set_ylabel('R [m]')
    axes[0].set_title('Shock Radius R(t)')

    # log-log
    axes[1].loglog(t_a * 1e3, R_a, 'b-', lw=2, label='Sedov-Taylor')
    axes[1].loglog(num_t * 1e3, num_R, 'ro', ms=6, label='Numerical')
    t_mid = t_a[len(t_a) // 2]
    R_mid = sedov.shock_radius(t_mid)
    axes[1].loglog(t_a * 1e3, R_mid * (t_a / t_mid) ** 0.4, 'k--',
                   alpha=0.4, label='∝ t^{2/5}')
    axes[1].set_xlabel('t [ms]')
    axes[1].set_ylabel('R [m]')
    axes[1].set_title('R(t) Log-Log')

    # Peak pressure
    axes[2].semilogy(t_a * 1e3, P_a / 1e6, 'b-', lw=2, label='Sedov-Taylor (R-H)')
    axes[2].semilogy(num_t * 1e3, num_P / 1e6, 'ro', ms=6, label='Numerical')
    axes[2].set_xlabel('t [ms]')
    axes[2].set_ylabel('Peak Pressure [MPa]')
    axes[2].set_title('Peak Pressure at Shock')

    for ax in axes:
        ax.legend()
        ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    path = os.path.join(output_dir, 'blast_shock_trajectory.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  {path}")


def plot_overpressure(params, output_dir):
    """過圧解析 (Rankine-Hugoniot, Hopkinson-Cranz)"""
    rh = RankineHugoniot(params.gamma)
    g = params.gamma
    p0 = params.p0

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle('Blast Overpressure Analysis', fontsize=14, fontweight='bold')

    # (1) Rankine-Hugoniot
    Ms = np.linspace(1.01, 10, 300)
    data = [rh.from_mach(m) for m in Ms]
    axes[0, 0].semilogy(Ms, [d['p_ratio'] for d in data], 'r-', lw=2, label='p₂/p₁')
    axes[0, 0].semilogy(Ms, [d['rho_ratio'] for d in data], 'b-', lw=2, label='ρ₂/ρ₁')
    axes[0, 0].semilogy(Ms, [d['T_ratio'] for d in data], 'g-', lw=2, label='T₂/T₁')
    axes[0, 0].set_xlabel('Shock Mach Number')
    axes[0, 0].set_ylabel('Ratio')
    axes[0, 0].set_title('Rankine-Hugoniot Relations')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3, which='both')

    # (2) 過圧 vs 換算距離
    Z = np.linspace(0.3, 50, 500)
    OP = p0 / 1e3 * (0.84 / Z + 2.7 / Z ** 2 + 7.05 / Z ** 3)
    axes[0, 1].semilogy(Z, OP, 'b-', lw=2)
    axes[0, 1].set_xlabel('Scaled Distance Z [m/kg^{1/3}]')
    axes[0, 1].set_ylabel('Peak Overpressure [kPa]')
    axes[0, 1].set_title('Kingery-Bulmash Scaling')
    for val, lbl, col in [(200, 'Heavy damage', 'darkred'),
                          (35, 'Structural', 'red'),
                          (7, 'Windows', 'orange'),
                          (2, 'Threshold', 'green')]:
        axes[0, 1].axhline(val, color=col, ls='--', alpha=0.5)
        axes[0, 1].text(42, val * 1.15, f'{lbl}\n({val} kPa)', fontsize=7, color=col)
    axes[0, 1].grid(True, alpha=0.3, which='both')

    # (3) 複数エネルギーの過圧 vs 距離
    W_list = [0.1, 1, 10, 100, 1000]
    cols = ['blue', 'green', 'orange', 'red', 'darkred']
    dist = np.linspace(1, 500, 1000)
    for W, col in zip(W_list, cols):
        Zd = dist / W ** (1 / 3)
        OP_d = p0 / 1e3 * (0.84 / Zd + 2.7 / Zd ** 2 + 7.05 / Zd ** 3)
        axes[1, 0].semilogy(dist, OP_d, '-', color=col, lw=1.5, label=f'{W} kg TNT')
    axes[1, 0].set_xlabel('Distance [m]')
    axes[1, 0].set_ylabel('Peak Overpressure [kPa]')
    axes[1, 0].set_title('Overpressure vs Distance')
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].set_xlim(1, 500)
    axes[1, 0].set_ylim(0.5, 1e4)
    axes[1, 0].grid(True, alpha=0.3, which='both')

    # (4) 等過圧線 (エネルギー vs 距離)
    W_g = np.logspace(-1, 4, 150)
    R_g = np.linspace(1, 500, 150)
    Wm, Rm = np.meshgrid(W_g, R_g)
    Zm = Rm / Wm ** (1 / 3)
    OPm = p0 / 1e3 * (0.84 / Zm + 2.7 / Zm ** 2 + 7.05 / Zm ** 3)
    lvls = [2, 7, 35, 100, 200, 500]
    cs = axes[1, 1].contour(Wm, Rm, OPm, levels=lvls, cmap='hot_r', linewidths=1.5)
    axes[1, 1].clabel(cs, fmt='%g kPa', fontsize=8)
    axes[1, 1].set_xlabel('TNT Equivalent [kg]')
    axes[1, 1].set_ylabel('Distance [m]')
    axes[1, 1].set_title('Overpressure Contours')
    axes[1, 1].set_xscale('log')
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(output_dir, 'blast_overpressure.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  {path}")


def plot_rt_contour(snapshots, params, output_dir):
    """r-t 圧力コンター図"""
    sedov = SedovTaylor(params)

    times = np.array([s['t'] for s in snapshots])
    r = snapshots[0]['r']
    r_lim = params.r_max * 0.5
    mask_r = r < r_lim
    r_plot = r[mask_r]
    P_data = np.array([s['p'][mask_r] for s in snapshots])

    fig, ax = plt.subplots(figsize=(10, 7))
    T_grid, R_grid = np.meshgrid(times * 1e3, r_plot, indexing='ij')
    cf = ax.contourf(R_grid, T_grid, np.log10(np.maximum(P_data, 1)),
                     levels=40, cmap='inferno')
    plt.colorbar(cf, ax=ax, label='log₁₀(p [Pa])')

    t_line = np.linspace(times[0], times[-1], 300)
    R_line = [sedov.shock_radius(t) for t in t_line]
    ax.plot(R_line, t_line * 1e3, 'w--', lw=2, label='Sedov-Taylor shock')

    ax.set_xlabel('r [m]')
    ax.set_ylabel('t [ms]')
    ax.set_title('Pressure Field r-t Diagram', fontweight='bold')
    ax.legend()
    ax.set_xlim(0, r_lim)
    plt.tight_layout()
    path = os.path.join(output_dir, 'blast_pressure_rt.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  {path}")


def print_shock_table(params):
    """解析テーブル出力"""
    sedov = SedovTaylor(params)
    rh = RankineHugoniot(params.gamma)
    W_kg = params.E0 / 4.184e6

    print("\n" + "=" * 88)
    print("衝撃波解析 (Sedov-Taylor + Rankine-Hugoniot)")
    print("=" * 88)
    print(f"E₀ = {params.E0:.3e} J ({W_kg:.1f} kg TNT),  "
          f"ρ₀ = {params.rho0} kg/m³,  p₀ = {params.p0:.0f} Pa,  "
          f"c₀ = {params.c0:.1f} m/s,  γ = {params.gamma}")
    print(f"Sedov α = {sedov.alpha:.4f}")
    print("-" * 88)
    print(f"{'t[ms]':>7} {'R[m]':>7} {'Vs[m/s]':>9} {'Ms':>6} "
          f"{'Δp[kPa]':>10} {'ρ₂/ρ₁':>7} {'T₂/T₁':>7}")
    print("-" * 88)

    for t_ms in [0.5, 1, 2, 5, 10, 15, 20, 30, 50, 80]:
        t = t_ms * 1e-3
        if t > params.t_end * 2:
            break
        Ms = sedov.shock_mach(t)
        if Ms < 1.001:
            continue
        R = sedov.shock_radius(t)
        Vs = sedov.shock_velocity(t)
        d = rh.from_mach(Ms)
        dp = (d['p_ratio'] - 1) * params.p0
        print(f"{t_ms:7.1f} {R:7.2f} {Vs:9.1f} {Ms:6.2f} "
              f"{dp/1e3:10.1f} {d['rho_ratio']:7.2f} {d['T_ratio']:7.2f}")

    print("=" * 88)


# =============================================================================
# メイン
# =============================================================================

def main():
    params = BlastParams(
        E0=4.184e9,     # 1 トン TNT
        rho0=1.225,
        p0=1.01325e5,
        gamma=1.4,
        r_max=150.0,
        Nr=3000,
        t_end=0.08,
        CFL=0.45,
    )

    output_dir = 'results'
    os.makedirs(output_dir, exist_ok=True)

    print_shock_table(params)

    print("\n" + "=" * 70)
    print("数値シミュレーション: 1D球対称 圧縮性オイラー方程式")
    print("手法: HLL Riemann Solver + Strang Splitting (幾何学的ソース項)")
    print("=" * 70)

    solver = EulerSolver(params)
    t_snap = np.array([1, 2, 5, 10, 15, 20, 30, 50, 80]) * 1e-3
    t_snap = t_snap[t_snap <= params.t_end]

    snapshots = solver.run(t_snap)

    print("\nグラフ生成中...")
    plot_profiles(snapshots, params, output_dir)
    plot_shock_trajectory(snapshots, params, output_dir)
    plot_overpressure(params, output_dir)
    if len(snapshots) >= 3:
        plot_rt_contour(snapshots, params, output_dir)

    print(f"\n結果 → '{output_dir}/':")
    for f in sorted(os.listdir(output_dir)):
        sz = os.path.getsize(os.path.join(output_dir, f)) / 1024
        print(f"  {f}  ({sz:.0f} KB)")


if __name__ == '__main__':
    main()
