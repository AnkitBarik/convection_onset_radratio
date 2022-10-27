#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
from scipy.special import airy
from scipy.io import loadmat

def get_airy_sol(s,chi,E):

    lamb = -1.02657 - 0.82534j
    rho0 = -2.338

    ri = chi/(1-chi)

    x = (s-ri)/E**(2./9.)

    Ai, Aip, Bi, Bip = airy(-lamb * x + rho0)

    Ai = np.abs(Ai)
    vs = Ai/max(Ai)

    return vs

mode = 'dark'

if mode == 'dark':
    plt.style.use("dark_background")
    plt.rcParams.update({
        "axes.facecolor"   : "#1b1b1b",
        "figure.facecolor" : "#1b1b1b",
        "figure.edgecolor" : "#1b1b1b",
        "savefig.facecolor": "#1b1b1b",
        "savefig.edgecolor": "#1b1b1b"})
    dormyCol = 'w'
else:
    dormyCol = 'k'

tksz = 30
axfontsize=40
lgfontsize=20
figW = 11
figH = 9
figsize = (figW, figH)
schematic=True

ekDirs = ['Ek1e-3','Ek1e-4','Ek1e-5',  'Ek1e-6','Ek1e-7','Ek3e-8','Ek1e-8','Ek3e-9','Ek1e-9']
ek = 10**np.array([-3,-4,-5,-6,-7,-7.5,-8,-8.5,-9])
chi = np.arange(0.05, 0.98, 0.03)
colors = cm.plasma_r(np.linspace(0.1, 0.8, len(ek)))
idx0_35 = np.argmin(np.abs(chi-0.35))

fig, ax = plt.subplots(1, 1, figsize=figsize)

chi1 = 0.35
chiStr = "{:.2f}".format(chi1)

for i in range(len(ekDirs)):
    os.chdir(ekDirs[i])
    os.chdir(chiStr)
    dat = np.loadtxt('vsEq.dat')
    s = dat[:, 0]/(1 - chi1)
    vs = dat[:, 1]
    ri = chi1/(1-chi1)
    if i == 0 or i == len(ek) - 1:
        label = r'$10^{%d}$' %(np.log10(ek[i]))
    else:
        label = None
    idx = np.argmin(np.abs(s-ri))

    vs_dormy = get_airy_sol(s,chi1,ek[i])

    ax.plot((s-ri)/ek[i]**(2./9.), vs/vs.max(), color=colors[i], label=label)
    ax.set_xlim(-1, 12)
    print(s[np.argmax(vs)])
    ax.set_ylim(-0.1, 1.1)
    print(os.getcwd())
    os.chdir('..')
    os.chdir('..')

label = 'Dormy et al., 2004'
ax.plot((s-ri)/ek[-1]**(2./9.), vs_dormy/vs_dormy.max(), color=dormyCol, label=label,ls='-',lw=5,alpha=0.6)
ax.set_xlabel(r'$(s-r_i)/E^{2/9}$', fontsize=1.5*axfontsize)
ax.set_ylabel(r'$|u_s|/|u_s|_{max}$', fontsize=1.5*axfontsize)
ax.tick_params(labelsize=tksz)
ax.legend(shadow=False,fontsize=1.5*lgfontsize,framealpha=0,
          frameon=True,edgecolor='none',
          loc='upper right',
          #title=r'$E$',title_fontsize=1.5*lgfontsize,
          borderpad=0.01)

plt.tight_layout()
plt.subplots_adjust(right=0.972)
# plt.show()
plt.savefig('dormySol_comp_'+chiStr+'.pdf',bbox_inches='tight')
