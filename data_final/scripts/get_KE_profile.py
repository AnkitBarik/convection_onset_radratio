#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from glob import glob
import os
from koreviz import *

def KE_sph_r(M):

    Cm = np.zeros([M.Q.shape[0], M.Q.shape[1]])
    norm = np.zeros_like(Cm)
    r2d, l2d = np.meshgrid(M.r, M.sh.l, indexing='ij')
    L = l2d * (l2d+1)
    norm = 4*np.pi/(2*l2d + 1)

    Cm[:, M.sh.m == 0] = 1.
    Cm[:, M.sh.m !=0]  = 2.

    KE_pol = Cm*norm*(np.abs(M.Q)**2 + L * np.abs(M.S)**2)
    KE_tor = Cm*norm * L * np.abs(M.T)**2
    KE_tot = KE_pol + KE_tor
    KE_pol_r = np.sum(KE_pol, axis=1)
    KE_tor_r = np.sum(KE_tor, axis=1)
    KE_r = KE_pol_r + KE_tor_r

    return KE_pol_r, KE_tor_r, KE_r

mdir = glob("Rac_m*")[0] #Replace with data directory

os.chdir(mdir)

pwd = os.getcwd()

sys.path.insert(0, pwd)

M = kmode(transform=False)
KE_pol_r, KE_tor_r, KE_r = KE_sph_r(M)

os.chdir('..')

dat = np.vstack([M.r, KE_pol_r, KE_tor_r, KE_r])

np.savetxt('Ekin_profile.dat', dat.transpose())
