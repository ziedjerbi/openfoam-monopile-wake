#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Estimation de l'energie cinetique turbulente totale dans le sillage du monopieu :
    k_total = <k_modelise> + k_resolu
    k_resolu = 1/2 (var Ux + var Uy + var Uz), variances sur les instants echantillonnes
Entrees : postProcessing/sampleDict/<t>/axeSillage_U.xy et axeSillage_k.xy
Sorties : figures/tke_totale.png, data/tke_totale_z1.txt
Usage   : python3 scripts/tke_totale.py
"""

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter

ROOT = Path(__file__).resolve().parent.parent
PP = ROOT / 'postProcessing' / 'sampleDict'
D, U0 = 0.5, 0.3162

# Lecture des instants
dirs = sorted([d for d in PP.iterdir() if d.is_dir()], key=lambda d: float(d.name))
times = [d.name for d in dirs]
x = np.loadtxt(dirs[0] / 'axeSillage_U.xy')[:, 0]
U = np.array([np.loadtxt(d / 'axeSillage_U.xy')[:, 1:4] for d in dirs])  # (nt, nx, 3)
K = np.array([np.loadtxt(d / 'axeSillage_k.xy')[:, 1] for d in dirs])    # (nt, nx)


def tke(idx):
    k_mod = K[idx].mean(axis=0)
    k_res = 0.5 * U[idx].var(axis=0, ddof=1).sum(axis=1)
    return k_mod, k_res, k_mod + k_res


xs = x / D
m = (xs >= 2.5) & (xs <= 14.0)


def fit(k):
    p = np.polyfit(np.log(xs[m]), np.log(k[m] / U0**2), 1)
    return p[0], np.exp(p[1])


def ecart4(expo, coef):
    return (0.6 * 4.0**(-1.1)) / (coef * 4.0**expo)


idx_all = np.arange(len(dirs))
k_mod, k_res, k_tot = tke(idx_all)

# Controle : <k_modelise> doit redonner data/k_moyen_z1.txt
ref = np.loadtxt(ROOT / 'data' / 'k_moyen_z1.txt')
ctrl_x = np.max(np.abs(ref[:, 0] - x))
ctrl_k = np.max(np.abs(ref[:, 1] - k_mod))

print(f"Instants utilises : {len(dirs)} ({times[0]} a {times[-1]} s)")
print(f"Controle k_moyen_z1.txt : ecart max x = {ctrl_x:.2e}, ecart max k = {ctrl_k:.2e}")
print()
print(f"{'Cas':<34}{'exposant':>10}{'coef':>9}{'ecart x*=4':>12}")
print(f"{'Rogan et al. (2016)':<34}{-1.100:>10.3f}{0.600:>9.3f}{'':>12}")
for nom, idx in [("k modelise (20 instants)", idx_all),
                 ("k total (20 instants)", idx_all),
                 ("k total (10 premiers)", idx_all[:10]),
                 ("k total (10 derniers)", idx_all[10:])]:
    km, kr, kt = tke(idx)
    k = km if nom.startswith("k modelise") else kt
    e, c = fit(k)
    print(f"{nom:<34}{e:>10.3f}{c:>9.4f}{ecart4(e, c):>11.2f}x")

i4 = np.argmin(np.abs(xs - 4.0))
print()
print(f"Part resolue a x*=4 : {100 * k_res[i4] / k_tot[i4]:.0f} % de k total")

# Donnees
np.savetxt(ROOT / 'data' / 'tke_totale_z1.txt', np.column_stack([x, k_mod, k_res, k_tot]),
           fmt='%.6f %.8f %.8f %.8f', header='x[m] k_modelise k_resolu k_total [m2/s2], z=-1.0 m, 20 instants')

# Figure
e_m, c_m = fit(k_mod)
e_t, c_t = fit(k_tot)
xr = np.linspace(2.5, 14.0, 200)

plt.rcParams['font.size'] = 11
fig, ax = plt.subplots(figsize=(9, 6))
ax.loglog(xr, 0.6 * xr**(-1.1), '-', lw=2.6, color='#111111', zorder=4,
          label=r'Rogan et al. (2016) : $0{,}6\;x^{*\,-1{,}10}$')
ax.loglog(xs[m], k_mod[m] / U0**2, 'o', ms=4, color='#0072B2', alpha=0.45, zorder=2,
          label='k modélisé (moyenne sur 20 instants)')
ax.loglog(xs[m], k_tot[m] / U0**2, 's', ms=4, color='#009E73', alpha=0.85, zorder=3,
          label='k total estimé = modélisé + résolu')
ax.loglog(xr, c_m * xr**e_m, '--', lw=1.6, color='#0072B2', zorder=5,
          label=f'Ajustement k modélisé : exposant {e_m:.2f}'.replace('.', ','))
ax.loglog(xr, c_t * xr**e_t, '--', lw=2.0, color='#D55E00', zorder=6,
          label=f'Ajustement k total : exposant {e_t:.2f}'.replace('.', ','))

ax.set_xlabel(r'Distance adimensionnée  $x^* = x/D$', fontsize=12)
ax.set_ylabel(r'$k\,/\,U_0^{2}$', fontsize=12)
ax.set_title("Sillage du monopieu — k modélisé et k total estimé (z = −1,0 m)\n"
             "Part résolue estimée sur 20 instants espacés de 5 s : estimation approchée",
             fontsize=11.5, pad=12)
ticks = [2.5, 3, 4, 5, 6, 8, 10, 14]
ax.xaxis.set_major_locator(FixedLocator(ticks))
ax.xaxis.set_major_formatter(FixedFormatter([str(t).replace('.', ',') for t in ticks]))
ax.xaxis.set_minor_locator(FixedLocator([]))
yt = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
ax.yaxis.set_major_locator(FixedLocator(yt))
ax.yaxis.set_major_formatter(FixedFormatter(['0,01', '0,02', '0,05', '0,10', '0,20', '0,50']))
ax.yaxis.set_minor_locator(FixedLocator([]))
ax.grid(True, which='major', alpha=0.25)
ax.legend(fontsize=9.5, loc='upper right', framealpha=0.95)
ax.set_xlim(2.4, 14.6)
ymax = max(0.3, 1.3 * np.max(k_tot[m] / U0**2))
ax.set_ylim(0.01, ymax)
fig.tight_layout()
fig.savefig(ROOT / 'figures' / 'tke_totale.png', dpi=220)
print(f"Figure : {ROOT / 'figures' / 'tke_totale.png'}")