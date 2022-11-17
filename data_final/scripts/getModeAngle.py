#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from glob import glob
import os
from koreviz import *

def getBlayerThicknessVal(s,KE,ri):
    sMax = s[np.argmax(KE)]
    mask = s >= sMax
    KE = KE[mask]
    s = s[mask]
    pcent = 1.

    s_end = s[np.argmin(np.abs(KE - KE.max()*pcent/100.))]

    return s_end - ri



def getBlayerThicknessD(r,dat,chi,ri,argrel):

    mask = r >= chi
    dat = dat[mask]
    r = r[mask]
    nr = len(r)
    dat= (dat/dat.max())
    if argrel:
        idx = argrelextrema(dat,np.greater)[0][1]
    else:
        idx = np.argmax(dat)
    A = np.array([ r[idx],dat[idx] ])
    B = np.array([ r[-1],dat[-1] ])
    b = B - A
    bhat = b/norm(b)
    norms = np.zeros(nr)

    for i in range(idx,nr):
        P = np.array([r[i],dat[i]])
        p = P - A
        norms[i] = norm( p - np.dot(p,b) * bhat)

    idx2 = np.argmax(norms)

    return (r[idx2] - ri)

dat = np.loadtxt('KE_zavg_s.dat')
s = dat[:,0]
KE_s = dat[:,1]

os.chdir(glob("Rac_m*")[0])

pwd = os.getcwd()

ek = pwd.split('/')[-2]

sys.path.insert(0, pwd)
M = kmode()

ri = M.r.min()
chi = ri
argrel = False
chiThick = ri

if ek == 1e-4 and chi in [0.05,0.08]:
    argrel = True
if ek == 1e-6 and chi == 0.05:
    chiThick = 0.19

radExt = getBlayerThicknessD(s,KE_s,chiThick,ri,argrel)

mask1 = M.r < (radExt + M.r.min())
half = M.ntheta//2
phi = M.phi[:M.nphi]
phiMax = np.zeros(M.nr)

for k in range(M.nr):
    phiMax[k] = phi[np.argmax(M.ur[k,half,:M.nphi])]

r = M.r[mask1]
phiM = phiMax[mask1]

idx1 = np.where(mask1)[0][0]
idx2 = np.argmax(phiM)

radExtFit = M.r[idx1 + idx2]
mask2 = M.r <= radExtFit

r = M.r[mask2]
phiM = phiMax[mask2]

x = r*np.cos(phiM)
y = r*np.sin(phiM)

X = np.vstack([x,y])

os.chdir('..')

np.savetxt('modeShape.dat',X.transpose())
