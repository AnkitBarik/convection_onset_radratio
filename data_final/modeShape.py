#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

plotStyle = 'light'

if plotStyle == 'dark':
    plt.style.use("dark_background")
    plt.rcParams.update({
        "axes.facecolor"   : "#1b1b1b",
        "figure.facecolor" : "#1b1b1b",
        "figure.edgecolor" : "#1b1b1b",
        "savefig.facecolor": "#1b1b1b",
        "savefig.edgecolor": "#1b1b1b"})
    spanCol = 'gray'
    schemCol = '#03fcb1'
    lineCol = 'w'
else:
    spanCol = 'gray'
    schemCol = '#03fcb1'
    lineCol = 'k'

tksz = 40
axfontsize=50
lgfontsize=30
figW = 11
figH = 9
figsize = (figW, figH)
figsizeSubplots = (27, 9)
figsizeE = (16, 9)

major_ticks = np.arange(0., 1.1, 0.2)
minor_ticks = np.arange(0.05, 0.96, 0.03)

def get_spirality(r,phi):

    delta_r = r.max() - r.min()
    delta_phi = phi.max() - phi.min()

    b_sp = np.abs(delta_phi)/np.abs(delta_r)

    #dy = np.gradient(y,x)
    #ddy = np.gradient(dy,x)

    #curvature =  np.abs(ddy)/( (1 + dy**2)**(3./2.) )

    return b_sp



def set_ax_params_nogrid(ax):
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xticks(major_ticks)
    ax.set_xlabel(r'$\chi = r_i/r_o$', fontsize=axfontsize)
    ax.tick_params(which='major', labelsize=tksz, length=10, direction='out')
    ax.tick_params(which='minor', length=5)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

ekDirs = ['Ek1e-4', 'Ek1e-5',  'Ek1e-6', 'Ek1e-7']
ek = [1e-4, 1e-5, 1e-6, 1e-7]
#ekDirs = ['Ek1e-6']
#ek = [1e-4]
chi = np.arange(0.05, 0.98, 0.03)

colors = cm.plasma_r(np.linspace(0.2, 0.8, 5))[1:]
fig, ax = plt.subplots(figsize=(14, 9))

b_sp = np.zeros_like(chi)

for i in range(len(ekDirs)):
    os.chdir(ekDirs[i])
    chiDirs = np.sort(next(os.walk('.'))[1])
    for j in range(len(chiDirs)):
        os.chdir(chiDirs[j])
        dat = np.loadtxt('modeShape.dat')
        
        r = dat[:,0]
        phi = dat[:,1]

        b_sp[j] = get_spirality(r,phi)

        if b_sp[j] <= 1e-4: # Some modes are pretty straight lines, so one gets near zero values
            b_sp[j] = np.nan

        print(os.getcwd())
        os.chdir('..')

    ax.semilogy(chi, b_sp, 'o', color=colors[i], label=r'$10^{%d}$' %np.log10(ek[i]))
    os.chdir('..')


set_ax_params_nogrid(ax)
ax.set_ylabel(r'$\frac{\Delta\phi}{\Delta s}$', fontsize=1.2*axfontsize, rotation=0, labelpad=50)
ax.legend(fontsize=lgfontsize, frameon=False)

plt.tight_layout()
#plt.show()
plt.savefig('modeSpiral.pdf',dpi=300)
