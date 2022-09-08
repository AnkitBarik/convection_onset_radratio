#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

tksz = 30
axfontsize=40
lgfontsize=20
figW = 11
figH = 9
figsize = (figW, figH)
schematic=True

ekDirs = ['Ek1e-3', 'Ek1e-4', 'Ek1e-5',  'Ek1e-6', 'Ek1e-7']
ek = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
chi = np.arange(0.05, 0.98, 0.03)
colors = cm.plasma_r(np.linspace(0.2, 0.8, 5))


if schematic:

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    for i in range(len(ekDirs)):
        os.chdir(ekDirs[i])
        os.chdir('0.50')
        dat = np.loadtxt('dissipation_profile.dat')
        r = dat[:, 0]
        Dnu = dat[:, 1]
        label = r'$10^{%d}$' %(np.log10(ek[i]))
        ax.plot((r-r[-1])/(ek[i])**(0.38), Dnu/Dnu[-1], color=colors[i], label=label)
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 1)
        print(os.getcwd())
        os.chdir('..')
        os.chdir('..')

    ax.set_xlabel(r'$(r-r_i)/E^{0.38}$', fontsize=axfontsize)
    ax.set_ylabel(r'$\mathcal{D}_{\nu}(r)/\mathcal{D}_{\nu}(r_i)$', fontsize=axfontsize)
    ax.tick_params(labelsize=tksz)
    ax.legend(fontsize=lgfontsize, frameon=False)

else:

    fig, axs = plt.subplots(6, 6, figsize=(12, 10))

    for i in range(len(ekDirs)):
        os.chdir(ekDirs[i])
        chiDirs = np.sort(next(os.walk('.'))[1])
        for j in range(31):
            os.chdir(str(chiDirs[j]))
            print(chiDirs[j])
            dat = np.loadtxt('dissipation_profile.dat')
            r = dat[:, 0]
            Dnu = dat[:, 1]
            jx = j//6
            jy = j%6
            axs[jx, jy].plot((r-r[-1])/(ek[i])**(0.38), Dnu/Dnu[-1], color=colors[i])
            axs[jx, jy].set_xlim(0, 5)
            axs[jx, jy].set_ylim(0, 0.3)
            print(os.getcwd())
            os.chdir('..')
        os.chdir('..')

plt.tight_layout()
plt.show()
#plt.savefig('Diss_icb.pdf',bbox_inches='tight')
