#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from numpy.linalg import norm
import os
from scipy.io import loadmat

plotStyle = 'light'

if plotStyle == 'dark':
    plt.style.use("dark_background")
    plt.rcParams.update({
        "axes.facecolor"   : "#1b1b1b",
        "figure.facecolor" : "#1b1b1b",
        "figure.edgecolor" : "#1b1b1b",
        "savefig.facecolor": "#1b1b1b",
        "savefig.edgecolor": "#1b1b1b"})
    spanCol = 'white'
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

    if idxMax == len(dat)-1:
        refVec = A
    elif idxMax == 0:
        refVec = B
    else:
        refVec = B

    for i in range(idx, nr):
        P = np.array([r[i], dat[i]])
        p = P - refVec
        norms[i] = norm( p - np.dot(p, b) * bhat)
    Idx = np.argmax(norms)

    return (r[Idx] - r.min())

def getBLayerMax(r, D):
    idx = np.argmax(D)
    rmin = r[idx]
    return rmin - r.min()

def getLineIntersect(m1, m2, c1, c2):

    x = (c2 - c1)/(m1 - m2)
    y = (m1*c2 - m2*c1)/(m1 - m2)

    return x, y

def getBlayerThicknessSlope(r,D,tolL=0.4,tolR=0.4,minD=False):

    if minD:
        d = getBLayerMax(r, D)
    else:
        d = getBlayerThicknessD(r, D)

    rBl = r.min() + d

    mask = (rFit < (rBl - tolL*d))

    p1 = np.polyfit(r[mask], D[mask], 1)

    mask = (rFit > (rBl + tolR*d))

    p2 = np.polyfit(r[mask], D[mask], 1)

    x, y = getLineIntersect(p1[0], p2[0], p1[1], p2[1])

    return x - r.min()

def plotSchematic(r,D,bound='icb'):

    D /= D.max()

    plt.figure(figsize=figsize)
    plt.plot(r, D, color=schemCol, lw=3, alpha=0.8)

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

    plt.axvspan(r.min(), x, alpha=0.3, color='gray')
    plt.ylim(-0.1, 1)
    plt.tick_params(labelsize=tksz)
    plt.xlabel(r'$r - r_i$', fontsize=axfontsize)
    plt.ylabel(r'$\mathcal{D}_{\nu}(r)/\mathcal{D}_{\nu}(r_i)$', fontsize=axfontsize)
    plt.tight_layout()
    plt.savefig('../../twoSlopeSchematic_'+ bound + plotStyle +'.pdf',dpi=300)
    # plt.close()

ekDirs = ['Ek1e-3', 'Ek1e-4', 'Ek1e-5',  'Ek1e-6', 'Ek1e-7']
ek = [1e-3,1e-4,1e-5,1e-6,1e-7]
chi = np.arange(0.05, 0.98, 0.03)

colors = cm.plasma_r(np.linspace(0.2, 0.8, 5))

delta_icb = np.zeros([len(ek), len(chi)])
ekExp = np.zeros(len(chi))
ekCoeff = np.zeros(len(chi))

for i in range(len(ekDirs)):
    os.chdir(ekDirs[i])
    chiDirs = np.sort(next(os.walk('.'))[1])
    for j in range(len(chiDirs)):
        os.chdir(str(chiDirs[j]))
        print(os.getcwd())
        dat = loadmat('uh_profiles.mat')
        r = dat['r'].flatten()/(1-chi[j])
        uh_mean = dat['uh_mean'].flatten()

        d = (r - r[-1])/(ek[i])**(0.4)
        mask = d <= 4.
        rFit  = (r[mask] - r[-1])
        uhFit = uh_mean[mask]

        if ekDirs[i] == 'Ek1e-6' and chiDirs[j] == '0.50':
            plotSchematic(rFit, uhFit, bound='icb')

        minD = False

        if ek[i] == 1e-3 and chi[j] == 0.14:
            "Changing tolerance!"
            tolL = 0.8
            tolR = 0.8
        else:
            tolL = 0.4
            tolR = 0.4
        tolL = 0.4
        tolR = 0.4

        delta_icb[i, j] = getBlayerThicknessSlope(rFit, uhFit, tolL=tolL, tolR=tolR, minD=minD)

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

delta_mean = np.zeros_like(chi)

def f(x, a):
    return x**a

fig, ax = plt.subplots(figsize=(14, 9))

for k in range(len(ek)):
    ax.semilogy(chi, delta_icb[k,:], '-o', color=colors[k], label=r'$10^{%d}$' %np.log10(ek[k]))

ax.set_ylabel(r'$\delta_{i}$', fontsize=axfontsize)
set_ax_params_nogrid(ax)
ax.legend(fontsize=lgfontsize, frameon=False, bbox_to_anchor=(1.01, 1.01), title=r'$E$', title_fontsize=lgfontsize)

plt.tight_layout()
plt.savefig('blthick_icb_' + plotStyle + '.pdf')

fig, ax = plt.subplots(figsize=figsize)

for k in range(len(chi)):
    p = np.polyfit(np.log10(ek[1:]), np.log10(delta_icb[1:, k]), 1)
    ekExp[k] = p[0]
    ekCoeff[k] = 10**p[1]

print("Scaling: a ek**b , a=%f, b=%f, std_dev=%f"  %(ekCoeff.mean(), ekExp.mean(), ekExp.std()))

ax.plot(chi, ekCoeff, '-o', label=r'$a$',color='#ff7f0e')
ax.plot(chi, ekExp, '-o', label=r'$b$',color='#1f77b4')
ax.axhline(y=1/3,lw=2,ls='--',color=lineCol,zorder=-9,label='Dormy et al., 2004')


set_ax_params_nogrid(ax)
ax.set_ylabel('Coefficients', fontsize=axfontsize)
ax.legend(fontsize=lgfontsize, frameon=False)
ax.set_yticks(np.arange(0, 3.2, 0.5))
ax.set_ylim(0, 2)

plt.tight_layout()
#plt.show()

plt.savefig('uhthick_icbCoeff_' + plotStyle + '.pdf',dpi=300)
