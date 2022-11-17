#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from numpy.linalg import norm
import os

plotStyle = 'light'
comp = 'full' # Can be full dissipation or toroidal ('tor')


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

def getBlayerThicknessD(r, dat):

    nr = len(r)
    dat= (dat/dat.max())
    idx = 0
    A = np.array([ r[idx], dat[idx] ])
    B = np.array([ r[-1], dat[-1] ])
    b = B - A
    bhat = b/norm(b)
    norms = np.zeros(nr)

    idxMax = np.argmax(dat)
    if idxMax == 0:
        refVec = A
    elif idxMax == len(dat)-1:
        refVec = B
    else:
        refVec = A

    for i in range(idx, nr):
        P = np.array([r[i], dat[i]])
        p = P - refVec
        norms[i] = norm( p - np.dot(p, b) * bhat)
    Idx = np.argmax(norms)

    return (r[Idx] - r.min())

def getLineIntersect(m1, m2, c1, c2):

    x = (c2 - c1)/(m1 - m2)
    y = (m1*c2 - m2*c1)/(m1 - m2)

    return x, y

def getBlayerThicknessSlope(r, D):

    d = getBlayerThicknessD(r, D)

    rBl = r.min() + d

    mask = (rFit < (rBl - 0.4*d))

    p1 = np.polyfit(r[mask], D[mask], 1)

    mask = (rFit > (rBl + 0.4*d))

    p2 = np.polyfit(r[mask], D[mask], 1)

    x, y = getLineIntersect(p1[0], p2[0], p1[1], p2[1])

    return x - r.min()

def plotSchematic(r,D):

    D /= D.max()

    plt.figure(figsize=figsize)
    plt.plot(r, D, color='#03fcb1', lw=3, alpha=0.8)

    d = getBlayerThicknessD(r, D)

    rBl = r.min() + d

    mask = (rFit < rBl - 0.4*d)

    p1 = np.polyfit(r[mask], D[mask], 1)

    plt.plot(r, np.polyval(p1, r), color=lineCol, lw=1.5, ls=':')

    mask = (rFit > rBl + 0.4*d)

    p2 = np.polyfit(r[mask], D[mask], 1)

    plt.plot(r, np.polyval(p2, r), color=lineCol, lw=1.5, ls=':')

    x, y = getLineIntersect(p1[0], p2[0], p1[1], p2[1])

    plt.plot(x, y, 'X',color=lineCol)

    plt.axvspan(r.min(), x, alpha=0.3, color=spanCol)

    plt.ylim(-0.1, 1)
    plt.tick_params(labelsize=tksz)
    plt.xlabel(r'$r_o - r$', fontsize=axfontsize)
    plt.ylabel(r'$\mathcal{D}_{\nu}(r)/\mathcal{D}_{\nu}(r_o)$', fontsize=axfontsize)
    plt.tight_layout()
    plt.savefig('../../twoSlopeSchematic_cmb' + plotStyle +'.pdf',dpi=300)
    # plt.close()

ekDirs = ['Ek1e-3', 'Ek1e-4', 'Ek1e-5',  'Ek1e-6', 'Ek1e-7']
ek = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
chi = np.arange(0.05, 0.98, 0.03)

colors = cm.plasma_r(np.linspace(0.2, 0.8, 5))

delta_cmb = np.zeros([len(ek), len(chi)])
ekExp = np.zeros(len(chi))
ekCoeff = np.zeros(len(chi))

for i in range(len(ekDirs)):
    os.chdir(ekDirs[i])
    chiDirs = np.sort(next(os.walk('.'))[1])
    for j in range(len(chiDirs)):
        os.chdir(str(chiDirs[j]))
        print(chiDirs[j])
        dat = np.loadtxt('dissipation_profile.dat')
        r = dat[:, 0]

        if comp == 'full':
            Dnu = dat[:, -1]/dat[0, -1]
        elif comp == 'tor':
            Dnu = dat[:, 2]/dat[0, 2]

        d = (r[0] - r)/np.sqrt(ek[i])
        mask = d <= 20
        rFit = r[0] - r[mask]
        DFit = Dnu[mask]
        if ekDirs[i] == 'Ek1e-6' and chiDirs[j] == '0.50':
            plotSchematic(rFit, DFit)
        delta_cmb[i, j] = getBlayerThicknessSlope(rFit, DFit)
        print(os.getcwd())
        os.chdir('..')

    os.chdir('..')

def set_ax_params_nogrid(ax):
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xticks(major_ticks)
    ax.set_xlabel(r'$\chi = r_i/r_o$', fontsize=axfontsize)
    ax.tick_params(which='major', labelsize=tksz, length=10, direction='out')
    ax.tick_params(which='minor', length=5)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


fig, ax = plt.subplots(figsize=(14, 9))

for k in range(len(ek)):
    ax.semilogy(chi, delta_cmb[k,:], '-o', color=colors[k], label=r'$10^{%d}$' %np.log10(ek[k]))

ax.set_ylabel(r'$\delta_{o}$', fontsize=axfontsize)
set_ax_params_nogrid(ax)
ax.legend(fontsize=lgfontsize, frameon=False, bbox_to_anchor=(1.01, 1.01), title=r'$E$', title_fontsize=lgfontsize)

plt.tight_layout()
plt.savefig('blthick_cmb_' + plotStyle + '.pdf')

fig, ax = plt.subplots(figsize=figsize)

for k in range(len(chi)):
    p = np.polyfit(np.log10(ek[1:]), np.log10(delta_cmb[1:, k]), 1)
    ekExp[k] = p[0]
    ekCoeff[k] = 10**p[1]

print("Scaling: a ek**b , a=%f, b=%f, std_dev=%f"  %(ekCoeff.mean(), ekExp.mean(), ekExp.std()))

ax.plot(chi, ekCoeff, '-o', label=r'$a$',color='#ff7f0e')
ax.plot(chi, ekExp, '-o', label=r'$b$',color='#1f77b4')
ax.axhline(y=0.5,lw=2,ls='--',color=lineCol,zorder=-9,label='Ekman layer')

set_ax_params_nogrid(ax)
ax.set_ylabel('Coefficients', fontsize=axfontsize)
ax.legend(fontsize=lgfontsize, frameon=False)
ax.set_yticks(np.arange(0, 3.2, 0.5))

plt.tight_layout()
#plt.show()

plt.savefig('blthick_cmbCoeff_' + plotStyle +  '.pdf',dpi=300)
