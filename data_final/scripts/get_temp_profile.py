#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from glob import glob
import os
import sys
from koreviz import *

def temp_sph_r(M):

    Cm   = np.zeros([M.Q.shape[0], M.Q.shape[1]])
    norm = np.zeros_like(Cm)
    r2d, l2d = np.meshgrid(M.r, M.sh.l, indexing='ij')

    L = l2d * (l2d+1)

    Cm[:, M.sh.m==0] = 1.0
    Cm[:, M.sh.m!=0] = 2.0

    norm = 4*np.pi/(2*l2d + 1)

    temp_l_r = Cm*norm*(np.abs(M.Q)**2)
    temp_r = np.sum(temp_l_r, axis=1)

    return np.sqrt(temp_r)

mdir = glob("Rac_m*")[0] #Replace with data directory

os.chdir(mdir)

pwd = os.getcwd()

sys.path.insert(0, pwd)

M = kmode(solnum=0, field='temp', transform=False)

temp_rms = temp_sph_r(M)

os.chdir('..')

dat = np.vstack([M.r, temp_rms])

np.savetxt('temp_rms_profile.dat', dat.transpose())
