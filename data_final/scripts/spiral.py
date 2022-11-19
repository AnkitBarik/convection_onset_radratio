#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from koreviz import *
import os
import sys

pwd = os.getcwd()

sys.path.insert(0,pwd)
sys.path.insert(0,pwd + '/bin')

import parameters as par

M = kmode(nthreads=4)
M.r /= (1-M.r.min())

smax = 5.*(par.Ek_gap)**(2/9) + M.r.min()
mask = M.r < smax
half = M.ntheta//2
ur = M.ur[mask,half,:]
r = M.r[mask]
ur_norm = np.zeros_like(ur)

for i in range(len(r)):
    ur_norm[i,:] = ur[i,:]/(np.abs(ur[i,:]).max())

phi = M.phi[:M.nphi]
phase = np.zeros(ur_norm.shape[0])

for i in range(len(r)):
    corr = np.correlate(ur_norm[i,:],ur_norm[0,:],'same')
    phase[i] = phi[np.argmax(corr)]

# Get rid of inner boundary
r = r[:-1]
phase = phase[:-1]

dPhi = np.diff(phase)
idx = np.argmax(np.abs(dPhi))

if np.abs(dPhi[idx]) > np.pi/M.m:
    if dPhi[idx] < 0:
        phi[idx+1:] += np.abs(dPhi[idx])
    else:
        phi[:idx] -= np.abs(dPhi[idx])

x = r*np.cos(phase)
y = r*np.sin(phase)

X = np.vstack([r,phase,x,y])

np.savetxt('modeShape.dat',X.transpose())
