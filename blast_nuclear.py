#!/usr/bin/env python3
"""
核爆発 爆風シミュレーション
============================
公開データに基づく核爆発の爆風解析

参考文献:
  - Glasstone & Dolan, "The Effects of Nuclear Weapons", 1977 (米国政府公開文書)
  - Taylor, G.I., "The Formation of a Blast Wave by a Very Intense Explosion", 1950
  - Sedov, L.I., "Similarity and Dimensional Methods in Mechanics", 1959
  - Brode, H.L., "Numerical Solutions of Spherical Blast Waves", 1955

歴史的核実験の公開データ:
  - Trinity (1945):   ~20 kT TNT
  - Little Boy (広島): ~15 kT TNT
  - Fat Man (長崎):   ~21 kT TNT
  - Tsar Bomba (1961): ~50 MT TNT (史上最大)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# =============================================================================
# 定数
# =============================================================================

KT_TO_JOULE = 4.184e12      # 1 kT TNT → J
MT_TO_JOULE = 4.184e15      # 1 MT TNT → J
ATM = 1.01325e5             # 標準気圧 [Pa]
RHO0 = 1.225                # 大気密度 [kg/m³]
GAMMA = 1.4                 # 比熱比
C0 = 340.3                  # 音速 [m/s]

# Sedov-Taylor 定数 (gamma=1.4, d=3)
ALPHA_ST = 0.851


# =============================================================================
# 核兵器パラメータ (公開データ)
# =============================================================================

WEAPONS = {
    'Little Boy': {
        'yield_kT': 15,
        'description': '広島 (1945/8/6), ウラン型',
        'burst_height_m': 600,   # 爆発高度 [m]
        'color': '#e74c3c',
    },
    'Fat Man': {
        'yield_kT': 21,
        'description': '長崎 (1945/8/9), プルトニウム型',
        'burst_height_m': 503,
        'color': '#e67e22',
    },
    'Trinity': {
        'yield_kT': 20,
        'description': '初の核実験 (1945/7/16), プルトニウム型',
        'burst_height_m': 0,     # 地表爆発 (塔上30mだが近似的に地表)
        'color': '#2ecc71',
    },
    'W-76': {
        'yield_kT': 100,
        'description': 'SLBM弾頭 (典型的戦略核)',
        'burst_height_m': 0,
        'color': '#3498db',
    },
    'B-83': {
        'yield_kT': 1200,
        'description': '米国最大の現役核爆弾',
        'burst_height_m': 0,
        'color': '#9b59b6',
    },
    'Tsar Bomba': {
        'yield_kT': 50000,
        'description': '史上最大 (1961/10/30), ソ連, 50 MT',
        'burst_height_m': 4000,
        'color': '#1a1a2e',
    },
}


# =============================================================================
# Sedov-Taylor 解析
# =============================================================================

def shock_radius(E_J, t, rho0=RHO0):
    """衝撃波半径 R(t) [m]"""
    return (E_J / (ALPHA_ST * rho0)) ** 0.2 * t ** 0.4

def shock_velocity(E_J, t, rho0=RHO0):
    """衝撃波速度 [m/s]"""
    return 0.4 * shock_radius(E_J, t, rho0) / t

def shock_mach(E_J, t):
    return shock_velocity(E_J, t) / C0

def shock_overpressure_RH(E_J, t):
    """Rankine-Hugoniot ピーク過圧 [Pa]"""
    Ms = shock_mach(E_J, t)
    if Ms < 1.0:
        return 0.0
    return ATM * (2 * GAMMA * Ms**2 - (GAMMA - 1)) / (GAMMA + 1) - ATM


# =============================================================================
# Glasstone-Dolan 過圧モデル (核爆発に特化した経験式)
# =============================================================================

def overpressure_glasstone(W_kT, R_m):
    """
    Glasstone-Dolan モデルによる核爆発ピーク過圧 [kPa]

    W_kT: 威力 [kT TNT]
    R_m:  爆心距離 [m]

    Hopkinson-Cranz スケーリング: Z = R / W_kg^(1/3) [m/kg^(1/3)]
    1 kT TNT = 10^6 kg TNT (質量換算)
    出典: "The Effects of Nuclear Weapons" (1977), Chapter III
    """
    W_kg = W_kT * 1e6  # kT → kg
    Z = R_m / (W_kg ** (1.0/3.0))  # 換算距離 [m/kg^(1/3)]
    Z = max(Z, 0.05)

    # Kinney-Graham (1985) / Brode 修正式
    # 有効範囲: Z > 0.05 m/kg^(1/3)
    op_kPa = ATM / 1e3 * (0.84 / Z + 2.7 / Z**2 + 7.05 / Z**3)
    return op_kPa


def overpressure_glasstone_array(W_kT, R_m):
    """配列対応版"""
    R_m = np.asarray(R_m, dtype=float)
    W_kg = W_kT * 1e6
    Z = R_m / (W_kg ** (1.0/3.0))
    Z = np.maximum(Z, 0.05)
    op_kPa = ATM / 1e3 * (0.84 / Z + 2.7 / Z**2 + 7.05 / Z**3)
    return op_kPa


# =============================================================================
# 核爆発特有の効果
# =============================================================================

def thermal_radius(W_kT, cal_per_cm2=10):
    """
    熱線到達距離 [m]

    cal_per_cm2: 熱線フルエンス [cal/cm²]
      10 cal/cm² → 第3度熱傷
      5 cal/cm² → 第2度熱傷
      3 cal/cm² → 第1度熱傷
    """
    # Glasstone-Dolan 近似: R_th ∝ W^0.41
    # 10 cal/cm² at 1 kT ≈ 500 m
    return 500 * (W_kT ** 0.41) / (cal_per_cm2 / 10) ** 0.5


def fireball_radius(W_kT):
    """火球最大半径 [m] (Glasstone-Dolan)"""
    # 空中爆発: R_fb = 34.4 * W^0.4  (W in kT)
    return 34.4 * W_kT ** 0.4


def crater_radius(W_kT):
    """クレーター半径 [m] (地表爆発の場合)"""
    # R_cr ≈ 20 * W^(1/3) (W in kT), 乾燥土壌
    return 20 * W_kT ** (1.0/3.0)


def damage_radii(W_kT):
    """各種被害半径 [m]"""
    R = np.linspace(10, 100000, 100000)
    OP = overpressure_glasstone_array(W_kT, R)

    results = {}
    thresholds = {
        'total_destruction':   140,     # 完全破壊 (20 psi)
        'severe_damage':       35,      # 重度損壊 (5 psi)
        'moderate_damage':     17,      # 中度損壊
        'light_damage':        7,       # 窓ガラス破損 (1 psi)
        'window_breakage':     3.5,     # 窓ガラスひび
    }

    for name, thresh_kPa in thresholds.items():
        idx = np.where(OP >= thresh_kPa)[0]
        if len(idx) > 0:
            results[name] = R[idx[-1]]
        else:
            results[name] = 0

    results['fireball'] = fireball_radius(W_kT)
    results['thermal_3rd'] = thermal_radius(W_kT, 10)
    results['thermal_2nd'] = thermal_radius(W_kT, 5)

    return results


# =============================================================================
# 1D球対称数値シミュレーション
# =============================================================================

class NuclearBlastSolver:
    """核爆発スケールの球対称爆風ソルバー"""

    def __init__(self, E_J, r_max, Nr=4000, CFL=0.4):
        self.E_J = E_J
        self.gamma = GAMMA
        self.r_max = r_max
        self.Nr = Nr
        self.CFL = CFL

        self.dr = r_max / Nr
        self.r = (np.arange(Nr) + 0.5) * self.dr

        # 初期条件
        rho = np.full(Nr, RHO0)
        u = np.zeros(Nr)
        p = np.full(Nr, ATM)

        r_init = max(10.0 * self.dr, fireball_radius(E_J / KT_TO_JOULE) * 0.1)
        vol_init = (4.0/3.0) * np.pi * r_init**3
        p_init = E_J * (GAMMA - 1) / vol_init
        mask = self.r <= r_init
        p[mask] = p_init

        self.U = np.array([rho, rho * u, p/(GAMMA-1) + 0.5*rho*u**2])

    def _prim(self, U):
        rho = np.maximum(U[0], 1e-20)
        u = U[1] / rho
        p = (GAMMA-1) * (U[2] - 0.5*rho*u**2)
        return rho, u, np.maximum(p, 1e-10)

    def _flux(self, rho, u, p):
        E = p/(GAMMA-1) + 0.5*rho*u**2
        return np.array([rho*u, rho*u**2+p, (E+p)*u])

    def _hll(self, rL, uL, pL, rR, uR, pR):
        cL = np.sqrt(GAMMA*pL/rL)
        cR = np.sqrt(GAMMA*pR/rR)
        SL = np.minimum(uL-cL, uR-cR)
        SR = np.maximum(uL+cL, uR+cR)
        FL = self._flux(rL, uL, pL)
        FR = self._flux(rR, uR, pR)
        UL = np.array([rL, rL*uL, pL/(GAMMA-1)+0.5*rL*uL**2])
        UR = np.array([rR, rR*uR, pR/(GAMMA-1)+0.5*rR*uR**2])
        d = SR - SL + 1e-30
        Fm = (SR*FL - SL*FR + SL*SR*(UR-UL)) / d
        return np.where(SL >= 0, FL, np.where(SR <= 0, FR, Fm))

    def _step(self, dt):
        rho, u, p = self._prim(self.U)
        rg = np.concatenate([[rho[0]], rho, [rho[-1]]])
        ug = np.concatenate([[-u[0]], u, [u[-1]]])
        pg = np.concatenate([[p[0]], p, [p[-1]]])
        Fe = self._hll(rg[:-1], ug[:-1], pg[:-1], rg[1:], ug[1:], pg[1:])
        self.U -= dt/self.dr * (Fe[:, 1:] - Fe[:, :-1])

        # 幾何学的ソース項
        rho, u, p = self._prim(self.U)
        r_safe = np.maximum(self.r, self.dr)
        inv_r = 2.0 / r_safe
        fac = 1.0 / (1.0 + inv_r * np.abs(u) * dt)
        self.U[0] -= dt * inv_r * rho * u * fac
        self.U[1] -= dt * inv_r * rho * u**2 * fac
        E_tot = self.U[2]
        self.U[2] -= dt * inv_r * (E_tot + p) * u * fac

        self.U[0] = np.maximum(self.U[0], 1e-20)
        rho, u, p = self._prim(self.U)
        self.U[2] = np.maximum(p, ATM*1e-6)/(GAMMA-1) + 0.5*rho*u**2

    def run(self, t_snapshots):
        t_snapshots = np.sort(t_snapshots)
        snaps = []
        t, n, si = 0.0, 0, 0
        while si < len(t_snapshots):
            rho, u, p = self._prim(self.U)
            c = np.sqrt(GAMMA * p / rho)
            dt = self.CFL * self.dr / np.max(np.abs(u) + c)
            if t + dt >= t_snapshots[si]:
                dt = max(t_snapshots[si] - t, 1e-15)
            self._step(dt)
            t += dt
            n += 1
            if t >= t_snapshots[si] - 1e-12:
                rho, u, p = self._prim(self.U)
                snaps.append({'t': t, 'r': self.r.copy(),
                              'rho': rho.copy(), 'u': u.copy(), 'p': p.copy()})
                above = p > ATM * 1.5
                rs = self.r[np.where(above)[0][-1]] if np.any(above) else 0
                print(f"  t={t:.4f} s | R_shock≈{rs:.0f} m | "
                      f"max(p)={np.max(p):.3e} Pa | n={n}")
                si += 1
        return snaps


# =============================================================================
# 可視化
# =============================================================================

def plot_overpressure_comparison(output_dir):
    """全兵器の過圧-距離比較"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    fig.suptitle('Nuclear Blast Overpressure Analysis\n'
                 '(Glasstone-Dolan Model, 1977)', fontsize=14, fontweight='bold')

    # (1) 過圧 vs 距離
    ax = axes[0, 0]
    R = np.linspace(100, 50000, 5000)
    for name, w in WEAPONS.items():
        OP = overpressure_glasstone_array(w['yield_kT'], R)
        ax.semilogy(R/1000, OP, '-', color=w['color'], lw=2,
                    label=f"{name} ({w['yield_kT']} kT)")
    ax.axhline(35, color='red', ls='--', alpha=0.4)
    ax.text(40, 38, 'Severe damage (35 kPa)', fontsize=7, color='red')
    ax.axhline(7, color='orange', ls='--', alpha=0.4)
    ax.text(40, 7.5, 'Light damage (7 kPa)', fontsize=7, color='orange')
    ax.set_xlabel('Distance [km]')
    ax.set_ylabel('Peak Overpressure [kPa]')
    ax.set_title('Overpressure vs Distance')
    ax.legend(fontsize=7, loc='upper right')
    ax.set_xlim(0.1, 100)
    ax.set_ylim(0.5, 1e6)
    ax.grid(True, alpha=0.3, which='both')

    # (2) 被害半径比較 (棒グラフ)
    ax = axes[0, 1]
    names = list(WEAPONS.keys())
    y_pos = np.arange(len(names))

    fb_r = [fireball_radius(WEAPONS[n]['yield_kT']) for n in names]
    total_r = [damage_radii(WEAPONS[n]['yield_kT'])['total_destruction']/1000 for n in names]
    severe_r = [damage_radii(WEAPONS[n]['yield_kT'])['severe_damage']/1000 for n in names]
    light_r = [damage_radii(WEAPONS[n]['yield_kT'])['light_damage']/1000 for n in names]
    thermal_r = [damage_radii(WEAPONS[n]['yield_kT'])['thermal_3rd']/1000 for n in names]

    bars = ax.barh(y_pos, light_r, color='#f39c12', alpha=0.4, label='Window breakage (7 kPa)')
    ax.barh(y_pos, severe_r, color='#e74c3c', alpha=0.5, label='Severe damage (35 kPa)')
    ax.barh(y_pos, total_r, color='#c0392b', alpha=0.7, label='Total destruction (140 kPa)')
    ax.barh(y_pos, [r/1000 for r in fb_r], color='white', alpha=0.9, label='Fireball')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel('Radius [km]')
    ax.set_title('Damage Radii Comparison')
    ax.legend(fontsize=7, loc='lower right')
    ax.grid(True, alpha=0.3, axis='x')

    # (3) 衝撃波半径の時間発展 (Sedov-Taylor)
    ax = axes[1, 0]
    for name, w in WEAPONS.items():
        E = w['yield_kT'] * KT_TO_JOULE
        # 強い衝撃波の有効時間
        t_max = 30.0  # 秒
        t = np.linspace(1e-4, t_max, 1000)
        R_t = np.array([shock_radius(E, ti) for ti in t])
        Ms_t = np.array([shock_mach(E, ti) for ti in t])
        # マッハ1以上のみ表示
        valid = Ms_t > 1.0
        if np.any(valid):
            ax.plot(t[valid], R_t[valid]/1000, '-', color=w['color'], lw=2,
                    label=f"{name} ({w['yield_kT']} kT)")
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Shock Radius [km]')
    ax.set_title('Shock Wave Propagation (Sedov-Taylor)')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # (4) 換算距離スケーリング (全兵器が1本の曲線に乗る)
    ax = axes[1, 1]
    Z = np.linspace(0.5, 200, 1000)
    OP = 6.895 * (6.7/Z**3 + 1.0/Z**1.5)  # Brode式

    ax.semilogy(Z, OP, 'k-', lw=3, label='Glasstone-Dolan\nUniversal Curve')

    # 各兵器のデータ点をプロット
    for name, w in WEAPONS.items():
        R_sample = np.array([1, 2, 5, 10, 20, 50]) * 1000  # m
        W_kg = w['yield_kT'] * 1e6
        Z_sample = R_sample / (W_kg ** (1/3))
        OP_sample = overpressure_glasstone_array(w['yield_kT'], R_sample)
        ax.semilogy(Z_sample, OP_sample, 'o', color=w['color'], ms=5,
                    label=f"{name}")

    ax.set_xlabel('Scaled Distance Z [m / kg^{1/3}]')
    ax.set_ylabel('Peak Overpressure [kPa]')
    ax.set_title('Hopkinson-Cranz Scaling\n(All yields collapse to one curve)')
    ax.legend(fontsize=7, ncol=2)
    ax.set_xlim(0.5, 200)
    ax.set_ylim(0.5, 1e5)
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    path = os.path.join(output_dir, 'nuclear_overpressure.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  {path}")


def plot_damage_map(output_dir):
    """Little Boy (広島) の被害半径マップ"""
    fig, ax = plt.subplots(figsize=(12, 12))

    W_kT = 15  # Little Boy
    dmg = damage_radii(W_kT)
    fb = fireball_radius(W_kT)

    # 同心円で被害範囲を表示
    zones = [
        (fb, 'Fireball', '#ffffff', 0.9),
        (dmg['total_destruction'], 'Total Destruction\n(140 kPa, 20 psi)', '#c0392b', 0.6),
        (dmg['severe_damage'], 'Severe Damage\n(35 kPa, 5 psi)', '#e74c3c', 0.4),
        (dmg['moderate_damage'], 'Moderate Damage\n(17 kPa, 2.5 psi)', '#e67e22', 0.3),
        (dmg['light_damage'], 'Light Damage / Windows\n(7 kPa, 1 psi)', '#f39c12', 0.2),
        (dmg['thermal_3rd'], '3rd Degree Burns\n(10 cal/cm²)', '#ff6b6b', 0.1),
    ]

    # 大きい順に描画
    zones_sorted = sorted(zones, key=lambda x: x[0], reverse=True)
    max_r = zones_sorted[0][0]

    for r_m, label, color, alpha in zones_sorted:
        circle = plt.Circle((0, 0), r_m/1000, color=color, alpha=alpha,
                           label=f'{label}: {r_m/1000:.1f} km')
        ax.add_patch(circle)
        if r_m > fb:
            ax.annotate(f'{r_m/1000:.1f} km',
                       xy=(r_m/1000*0.707, r_m/1000*0.707),
                       fontsize=8, ha='center',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    ax.set_xlim(-max_r/1000*1.15, max_r/1000*1.15)
    ax.set_ylim(-max_r/1000*1.15, max_r/1000*1.15)
    ax.set_aspect('equal')
    ax.set_xlabel('Distance [km]')
    ax.set_ylabel('Distance [km]')
    ax.set_title(f'Little Boy (Hiroshima) — {W_kT} kT TNT\nBlast & Thermal Damage Radii',
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.2)

    # 爆心マーク
    ax.plot(0, 0, 'r*', ms=15, zorder=10)
    ax.text(0.1, 0.1, 'Ground Zero', fontsize=9, fontweight='bold')

    plt.tight_layout()
    path = os.path.join(output_dir, 'nuclear_damage_map_hiroshima.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  {path}")


def plot_numerical_profiles(output_dir):
    """Little Boy の数値シミュレーション結果"""
    W_kT = 15
    E_J = W_kT * KT_TO_JOULE

    print(f"\n数値シミュレーション: Little Boy ({W_kT} kT)")
    print(f"  E = {E_J:.3e} J")

    solver = NuclearBlastSolver(
        E_J=E_J, r_max=10000, Nr=5000, CFL=0.4)

    t_snap = np.array([0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0])
    snapshots = solver.run(t_snap)

    # プロット
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle(f'Nuclear Blast Wave — Little Boy ({W_kT} kT TNT)\n'
                 f'1D Spherical Euler Equations + HLL Solver',
                 fontsize=13, fontweight='bold')

    colors = plt.cm.inferno(np.linspace(0.15, 0.85, len(snapshots)))

    for i, s in enumerate(snapshots):
        c = colors[i]
        lbl = f't={s["t"]:.3f} s' if s['t'] < 1 else f't={s["t"]:.1f} s'
        axes[0].plot(s['r']/1000, s['rho'], '-', color=c, lw=1.2, label=lbl)
        axes[1].plot(s['r']/1000, s['u'], '-', color=c, lw=1.2, label=lbl)
        axes[2].plot(s['r']/1000, s['p']/1e6, '-', color=c, lw=1.2, label=lbl)

    for ax, yl, tt in zip(axes,
        ['Density [kg/m³]', 'Velocity [m/s]', 'Pressure [MPa]'],
        ['ρ(r)', 'u(r)', 'p(r)']):
        ax.set_xlabel('Distance [km]')
        ax.set_ylabel(yl)
        ax.set_title(tt)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 8)

    plt.tight_layout()
    path = os.path.join(output_dir, 'nuclear_profiles.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  {path}")

    return snapshots


def print_damage_table():
    """全兵器の被害半径テーブル"""
    print("\n" + "=" * 100)
    print("核爆発 被害半径一覧 (Glasstone-Dolan Model)")
    print("=" * 100)
    print(f"{'兵器':>15} {'威力':>10} {'火球':>8} {'完全破壊':>8} "
          f"{'重度損壊':>8} {'軽度損壊':>8} {'3度熱傷':>8} {'2度熱傷':>8}")
    print(f"{'':>15} {'[kT]':>10} {'[m]':>8} {'[km]':>8} "
          f"{'[km]':>8} {'[km]':>8} {'[km]':>8} {'[km]':>8}")
    print("-" * 100)

    for name, w in WEAPONS.items():
        kT = w['yield_kT']
        fb = fireball_radius(kT)
        d = damage_radii(kT)
        print(f"{name:>15} {kT:>10,} {fb:>8.0f} "
              f"{d['total_destruction']/1000:>8.1f} "
              f"{d['severe_damage']/1000:>8.1f} "
              f"{d['light_damage']/1000:>8.1f} "
              f"{d['thermal_3rd']/1000:>8.1f} "
              f"{d['thermal_2nd']/1000:>8.1f}")

    print("=" * 100)
    print("出典: Glasstone & Dolan, 'The Effects of Nuclear Weapons', 1977 (米国政府公開文書)")


def print_sedov_table():
    """衝撃波パラメータ (Little Boy)"""
    W_kT = 15
    E_J = W_kT * KT_TO_JOULE

    print(f"\n{'='*90}")
    print(f"Sedov-Taylor 衝撃波解析: Little Boy ({W_kT} kT)")
    print(f"{'='*90}")
    print(f"{'t [s]':>8} {'R [km]':>8} {'Vs [m/s]':>10} {'Ms':>7} "
          f"{'Δp [kPa]':>12} {'ρ₂/ρ₁':>7} {'T₂/T₁':>8}")
    print("-" * 90)

    times = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    for t in times:
        R = shock_radius(E_J, t)
        Vs = shock_velocity(E_J, t)
        Ms = Vs / C0
        if Ms < 1.001:
            continue
        g = GAMMA
        p_ratio = (2*g*Ms**2 - (g-1)) / (g+1)
        rho_ratio = ((g+1)*Ms**2) / ((g-1)*Ms**2 + 2)
        T_ratio = p_ratio / rho_ratio
        dp = (p_ratio - 1) * ATM
        print(f"{t:8.3f} {R/1000:8.2f} {Vs:10.0f} {Ms:7.1f} "
              f"{dp/1e3:12.0f} {rho_ratio:7.2f} {T_ratio:8.1f}")

    print(f"{'='*90}")


# =============================================================================
# VTK エクスポート (ParaView用)
# =============================================================================

def export_vtk(snapshots, output_dir):
    """数値シミュレーション結果をVTKエクスポート"""
    try:
        import vtk
        from vtk.util.numpy_support import numpy_to_vtk
    except ImportError:
        print("  VTK未インストール — スキップ")
        return

    vtk_dir = os.path.join(output_dir, 'vtk_nuclear')
    os.makedirs(vtk_dir, exist_ok=True)

    N = 80
    times = []
    files = []

    for i, s in enumerate(snapshots):
        r_1d, p_1d, rho_1d, u_1d = s['r'], s['p'], s['rho'], s['u']
        L = min(s['r'][-1] * 0.8, 8000)
        Nz = N // 2

        x = np.linspace(-L, L, N)
        y = np.linspace(-L, L, N)
        z = np.linspace(0, L, Nz)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        R = np.sqrt(X**2 + Y**2 + Z**2)

        p_3d = np.interp(R.ravel(), r_1d, p_1d).reshape(R.shape)
        rho_3d = np.interp(R.ravel(), r_1d, rho_1d).reshape(R.shape)

        rgrid = vtk.vtkRectilinearGrid()
        rgrid.SetDimensions(N, N, Nz)

        for arr_np, setter in [(x, rgrid.SetXCoordinates),
                                (y, rgrid.SetYCoordinates),
                                (z, rgrid.SetZCoordinates)]:
            va = vtk.vtkFloatArray()
            for v in arr_np:
                va.InsertNextValue(v)
            setter(va)

        for vtk_name, data_3d in [('Pressure', p_3d),
                                   ('Overpressure', p_3d - ATM),
                                   ('Density', rho_3d),
                                   ('LogPressure', np.log10(np.maximum(p_3d, 1)))]:
            arr = numpy_to_vtk(data_3d.ravel(order='F'), deep=True)
            arr.SetName(vtk_name)
            rgrid.GetPointData().AddArray(arr)

        rgrid.GetPointData().SetActiveScalars('Pressure')

        fpath = os.path.join(vtk_dir, f'nuclear_{i:04d}.vtr')
        writer = vtk.vtkXMLRectilinearGridWriter()
        writer.SetFileName(fpath)
        writer.SetInputData(rgrid)
        writer.SetDataModeToBinary()
        writer.Write()

        times.append(s['t'])
        files.append(fpath)
        print(f"  VTK [{i+1}/{len(snapshots)}]: t={s['t']:.3f} s")

    # PVD
    pvd_path = os.path.join(output_dir, 'nuclear_blast.pvd')
    lines = ['<?xml version="1.0"?>',
             '<VTKFile type="Collection" version="0.1">',
             '  <Collection>']
    for t, f in zip(times, files):
        lines.append(f'    <DataSet timestep="{t:.6e}" file="vtk_nuclear/{os.path.basename(f)}"/>')
    lines += ['  </Collection>', '</VTKFile>']
    with open(pvd_path, 'w') as f:
        f.write('\n'.join(lines))
    print(f"  PVD: {pvd_path}")


# =============================================================================
# メイン
# =============================================================================

def main():
    output_dir = 'results'
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("核爆発 爆風シミュレーション")
    print("=" * 60)

    # テーブル出力
    print_damage_table()
    print_sedov_table()

    # 可視化
    print("\nグラフ生成中...")
    plot_overpressure_comparison(output_dir)
    plot_damage_map(output_dir)

    # 数値シミュレーション (Little Boy)
    snapshots = plot_numerical_profiles(output_dir)

    # VTK エクスポート (ParaView用)
    print("\nParaView用 VTK エクスポート...")
    export_vtk(snapshots, output_dir)

    print(f"\n結果 → '{output_dir}/':")
    for f in sorted(os.listdir(output_dir)):
        if f.startswith('nuclear'):
            sz = os.path.getsize(os.path.join(output_dir, f)) / 1024
            print(f"  {f}  ({sz:.0f} KB)")

    print(f"\nParaView: paraview results/nuclear_blast.pvd")


if __name__ == '__main__':
    main()
