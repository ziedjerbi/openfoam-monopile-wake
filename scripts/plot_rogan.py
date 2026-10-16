#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comparaison de l'energie cinetique turbulente (TKE) simulee dans le sillage
d'un monopieu (URANS k-omega SST, OpenFOAM v2112) avec la loi experimentale
de Rogan et al. (2016), Renewable Energy 93:180-187.

Entree : data/k_moyen_z1.txt   (colonnes : x [m], k [m2/s2])
Sortie : figures/comparaison_rogan.png
Usage  : python3 scripts/plot_rogan.py   (depuis n'importe quel dossier)
"""

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter

# Dossier racine du depot (parent du dossier scripts/)
ROOT = Path(__file__).resolve().parent.parent
FICHIER_K = ROOT / 'data' / 'k_moyen_z1.txt'
FICHIER_FIG = ROOT / 'figures' / 'comparaison_rogan.png'

D    = 0.5        # diametre du monopieu, echelle modele (m)
U0   = 0.3162     # vitesse incidente (m/s)
K_IN = 3.75e-4    # k impose en entree (m2/s2), voir case/0/k

x, k = np.loadtxt(FICHIER_K, unpack=True)
xs  = x / D
tke = k / U0**2

m = (xs >= 2.5) & (xs <= 14.0)
xs, tke = xs[m], tke[m]

# Loi experimentale de Rogan et al. (2016)
xs_ref  = np.linspace(2.5, 14.0, 200)
tke_ref = 0.6 * xs_ref**(-1.1)

# Ajustement en loi de puissance sur la simulation
p = np.polyfit(np.log(xs), np.log(tke), 1)
expo, coef = p[0], np.exp(p[1])
tke_fit = coef * xs_ref**expo

# Niveau de turbulence impose en entree de la simulation
tke_in = K_IN / U0**2

plt.rcParams['font.size'] = 11
fig, ax = plt.subplots(figsize=(9, 6))

lab_fit = (rf'Ajustement : ${coef:.3f}\;x^{{*\,{expo:.2f}}}$'
           .replace('0.', '0{,}'))
lab_in = (f"Entrée de la simulation : $k/U_0^2$ = {tke_in:.4f}"
          .replace('.', ','))

ax.loglog(xs_ref, tke_ref, '-', lw=2.6, color='#111111', zorder=4,
          label=r'Rogan et al. (2016) : $0{,}6\;x^{*\,-1{,}10}$   ($R^2$ = 0,945)')
ax.loglog(xs, tke, 'o', ms=4.5, color='#0072B2', alpha=0.8, zorder=3,
          label='URANS $k$-$\\omega$ SST \u2014 20 instants moyennés')
ax.loglog(xs_ref, tke_fit, '--', lw=2.0, color='#D55E00', zorder=5,
          label=lab_fit)
ax.axhline(tke_in, ls=':', lw=1.6, color='#777777', zorder=2,
           label=lab_in)

# Fleche d'ecart a x* = 4
xa = 4.0
ya1, ya2 = coef * xa**expo, 0.6 * xa**(-1.1)
ax.annotate('', xy=(xa, ya2), xytext=(xa, ya1),
            arrowprops=dict(arrowstyle='<->', color='#444444', lw=1.4))
ax.text(xa * 1.06, np.sqrt(ya1 * ya2), f'\u00d7{ya2/ya1:.1f}'.replace('.', ','),
        fontsize=10, color='#444444', va='center', fontweight='bold')

ax.set_xlabel(r'Distance adimensionnée  $x^* = x/D$', fontsize=12)
ax.set_ylabel(r'$k\,/\,U_0^{2}$', fontsize=12)
ax.set_title("Sillage turbulent d'un monopieu offshore \u2014 comparaison URANS / Rogan et al. (2016)\n"
             r"Similitude de Froude 1:10   ·   $D$ = 0,5 m   ·   $U_0$ = 0,316 m/s"
             "   ·   643 000 cellules",
             fontsize=11.5, pad=12)

ticks = [2.5, 3, 4, 5, 6, 8, 10, 14]
labs  = ['2,5', '3', '4', '5', '6', '8', '10', '14']
ax.xaxis.set_major_locator(FixedLocator(ticks))
ax.xaxis.set_major_formatter(FixedFormatter(labs))
ax.xaxis.set_minor_locator(FixedLocator([]))

yt = [0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2]
ax.yaxis.set_major_locator(FixedLocator(yt))
ax.yaxis.set_major_formatter(FixedFormatter(['0,005', '0,01', '0,02', '0,03', '0,05', '0,10', '0,20']))
ax.yaxis.set_minor_locator(FixedLocator([]))

ax.grid(True, which='major', alpha=0.25, ls='-')
ax.legend(fontsize=9.5, loc='upper right', framealpha=0.95)
ax.set_xlim(2.4, 14.6)
ax.set_ylim(0.003, 0.30)

fig.tight_layout()
fig.savefig(FICHIER_FIG, dpi=220)
print(f"Figure enregistree : {FICHIER_FIG}")
print(f"Exposant ajuste : {expo:.3f}   (Rogan : -1.100)")
print(f"Coefficient     : {coef:.4f}   (Rogan :  0.600)")
print(f"Ecart a x*=4    : facteur {ya2/ya1:.2f}")
print(f"Niveau d'entree : k/U0^2 = {tke_in:.4f}")