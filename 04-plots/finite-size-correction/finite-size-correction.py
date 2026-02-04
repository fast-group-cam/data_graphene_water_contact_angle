#! /usr/bin/env python

import os
import configparser
import numpy as np
import matplotlib.pyplot as plt

PATTERN = '../../01-unstrained/mace/$$$$-molecules/run_prod/contact-angle/results.ini'
NUMBERS = [4680, 4300, 4000, 3500, 3000, 2500, 2000, 1500, 1000, 800, 600, 500, 300]
RAD_MIN = 11.9
RAD_MAX = 84.9

DATA_MA = [[0.012804348602725568, -0.0052326929],
           [0.027536317593098947, -0.0178014116],
           [0.03206243146910992, -0.0257763506],
           [0.043356454064036354, -0.0518425524]]

DATA_SPCE = [[0.0225372095, 0.1275552204],
             [0.025098998, 0.1114790951],
             [0.0301432591, 0.1027633133],
             [0.0377419546, 0.098671345]]

def cos_deg(t):
    return np.cos(t * np.pi / 180.0)

def cos_deg_err(t, err):
    return err * np.abs(np.sin(t * np.pi / 180.0)) * np.pi / 180.0

def arccos_deg(c):
    return np.arccos(np.clip(c, a_min=-1.0, a_max=1.0)) * 180.0 / np.pi

N_water = list()
rad = list()
rad_err = list()
ang = list()
ang_err = list()

for number in NUMBERS:
    filename = PATTERN.replace('$$$$', str(number))
    if os.path.isfile(filename):
        config = configparser.ConfigParser()
        config.read(filename)
        if 'Block-Averaged Interface' in config:
            a = config.getfloat('Block-Averaged Interface', 'Three-phase line radius [A], mean of block means', fallback=None)
            b = config.getfloat('Block-Averaged Interface', 'Three-phase line radius [A], uncertainty', fallback=None)
            c = config.getfloat('Block-Averaged Interface', 'Contact angle [deg], mean of block means', fallback=None)
            d = config.getfloat('Block-Averaged Interface', 'Contact angle [deg], uncertainty', fallback=None)
            if None not in (a, b, c, d):
                N_water.append(int(number))
                rad.append(a)
                rad_err.append(b)
                ang.append(c)
                ang_err.append(d)
    
N_water = np.array(N_water)
rad = np.array(rad)
rad_err = np.array(rad_err)
ang = np.array(ang)
ang_err = np.array(ang_err)
x = 1.0 / rad
y = cos_deg(ang)
y_err = cos_deg_err(ang, ang_err)

p, cov = np.polyfit(x, y, 1, w=np.power(y_err, -1), cov=True)
#p, cov = np.polyfit(x, y, 1, cov=True)
theta_inf = arccos_deg(p[1])
theta_inf_err = 2 * np.sqrt(cov[1,1] / (1.0 - (p[1]**2))) * 180 / np.pi
print(f'Best-fit gradient: {p[0]} \u00b1 {np.sqrt(cov[0,0])}')

hor_range = np.linspace(RAD_MIN, RAD_MAX, 200)

data_ma = np.array(DATA_MA)
data_spce = np.array(DATA_SPCE)

fig, ax = plt.subplots()
fig.set_size_inches(6.4, 4)

p_ma, cov_ma = np.polyfit(data_ma[:,0], data_ma[:,1], 1, cov=True)
theta_inf_ma = arccos_deg(p_ma[1])
theta_inf_err_ma = np.sqrt(cov_ma[1,1] / (1.0 - (p_ma[1]**2))) * 180 / np.pi
ax.plot(hor_range, arccos_deg(p_ma[1] + (p_ma[0] / hor_range)), '--', color='forestgreen')
ax.plot(1.0 / data_ma[:,0], arccos_deg(data_ma[:,1]), 's', color='forestgreen', label='Ma et al., TIP4P + LJ')

p_spce, cov_spce = np.polyfit(data_spce[:,0], data_spce[:,1], 1, cov=True)
theta_inf_spce = arccos_deg(p_spce[1])
theta_inf_err_spce = np.sqrt(cov_spce[1,1] / (1.0 - (p_spce[1]**2))) * 180 / np.pi
ax.plot(hor_range, arccos_deg(p_spce[1] + (p_spce[0] / hor_range)), '--', color='salmon')
ax.plot(1.0 / data_spce[:,0], arccos_deg(data_spce[:,1]), 'D', color='salmon', label='Carlson et al., SPC/E + LJ')

ax.plot(hor_range, arccos_deg(p[1] + (p[0] / hor_range)), '--', color='steelblue')
ax.errorbar(rad, ang, yerr=(2 * ang_err), xerr=(2 * rad_err), color='steelblue', fmt='o', mfc='none', mec='steelblue', label='This work, MLP (revPBE-D3)')
ax.set_xlabel(r'Three-phase contact line radius $a\;[\AA]$')
ax.set_ylabel(r'Microscopic contact angle $\theta\;[\degree]$')
ax.set_xlim(RAD_MIN, RAD_MAX)
ax.set_ylim(65.1, 102.9)

grey = (0.5, 0.5, 0.5)
N_mol = np.max(N_water)
N_at = (3 * N_mol) + 8640
ax.annotate(f'{N_mol:,d} molecules\n({N_at:,d} atoms)', (rad[0], ang[0]), (49.7, 77.9), arrowprops={'arrowstyle': '->', 'color': grey, 'shrinkB': 10}, color=grey, ha='left', va='center')
N_mol = np.min(N_water)
ax.annotate(f'{N_mol:,d} molecules', (rad[-1], ang[-1]), (13.5, 66.9), arrowprops={'arrowstyle': '->', 'color': grey, 'shrinkB': 22}, color=grey, ha='left', va='center')
ax.text(67.9, 89.3, r'$\theta_{\infty} = ' + f'{theta_inf_ma:.1f}' + r' \pm ' + f'{theta_inf_err_ma:.1f}' + r'\degree$', color='forestgreen', ha='left', va='top')
ax.text(67.9, 81.4, r'$\theta_{\infty} = ' + f'{theta_inf_spce:.1f}' + r' \pm ' + f'{theta_inf_err_spce:.1f}' + r'\degree$', color='salmon', ha='left', va='top')
ax.text(67.9, 72.3, r'$\theta_{\infty} = ' + f'{theta_inf:.1f}' + r' \pm ' + f'{theta_inf_err:.1f}' + r'\degree$', color='steelblue', ha='left', va='top')
ax.legend()

fig.savefig('finite-size.png', dpi=(4.558*fig.dpi), bbox_inches='tight', pad_inches=0.05)
