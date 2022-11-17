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

ekDirs = ['Ek1e-4', 'Ek1e-5',  'Ek1e-6', 'Ek1e-7']
ek = [1e-4, 1e-5, 1e-6, 1e-7]
chi = np.arange(0.05, 0.98, 0.03)
colors = cm.plasma_r(np.linspace(0.2, 0.8, 5)[1:])

if schematic:

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    chi1 = 0.92
    chiStr = "{:.2f}".format(chi1)

    for i in range(len(ekDirs)):
        os.chdir(ekDirs[i])
        os.chdir(chiStr)
        dat = np.loadtxt('KE_zavg_s.dat')
        s = dat[:, 0]/(1 - chi1)
        KE = dat[:, 1]
        ri = chi1/(1-chi1)
        label = r'$10^{%d}$' %(np.log10(ek[i]))
        idx = np.argmin(np.abs(s-ri))
        ax.plot((s-ri)/ek[i]**(2./9.), KE/KE.max(), color=colors[i], label=label)
        ax.set_xlim(0, 15)
        print(s[np.argmax(KE)])
        ax.set_ylim(-0.1, 1)
        print(os.getcwd())
        os.chdir('..')
        os.chdir('..')

    ax.set_xlabel(r'$(s-r_i)/E^{2/9}$', fontsize=1.5*axfontsize)
    ax.set_ylabel(r'$\mathcal{E}_{kin}(s)/\mathcal{E}_{kin_{max}}$', fontsize=1.5*axfontsize)
    ax.tick_params(labelsize=tksz)
    if chi1 == 0.11:
        ax.legend(shadow=False,fontsize=1.5*lgfontsize,framealpha=0,
                 frameon=True,edgecolor='none',
                 loc='upper right',
                 title=r'$E$',title_fontsize=1.5*lgfontsize,
                 borderpad=0.01)

else:

    fig, axs = plt.subplots(6, 6, figsize=(12, 10))

    for i in range(len(ekDirs)):
        os.chdir(ekDirs[i])
        chiDirs = np.sort(next(os.walk('.'))[1])
        for j in range(len(chiDirs)):
            os.chdir(str(chiDirs[j]))
            print(chiDirs[j])
            dat = np.loadtxt('KE_zavg_s.dat')
            s = dat[:, 0]/(1-chi[j])
            ri = chi[j]/(1-chi[j])
            KE = dat[:, 1]
            jx = j//6
            jy = j%6
            ekfac = ek[i]**(2/9)
            axs[jx, jy].plot((s-ri)/ekfac, KE/KE.max(), color=colors[i])
            axs[jx, jy].set_xlim(-1, 15)
            axs[jx, jy].set_ylim(0, 1)
            axs[jx, jy].set_title('%.2f' %chi[j])
            # s_end = getBlayerThicknessVal(s,KE)
            # axs[jx,jy].axvline((s_end-ri)/ekfac)
            print(os.getcwd())
            os.chdir('..')
        os.chdir('..')

plt.tight_layout()
# plt.show()
plt.savefig('radExtent_example_'+chiStr+'.pdf',bbox_inches='tight')
