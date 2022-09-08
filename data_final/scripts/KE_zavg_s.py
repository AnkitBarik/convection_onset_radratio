#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

from koreviz import *
import numpy as np
import os
from magic import *
from glob import glob

mdir = glob("Rac_m*")[0] #Replace with data directory

os.chdir(mdir)

pwd = os.getcwd()

sys.path.insert(0, pwd)

M    = kmode(field='u', nthreads=20)
half = int(M.ntheta/2)
KE   = 0.5*np.sqrt(M.ur**2 + M.utheta**2 + M.uphi**2)
KE2d = KE.mean(axis=-1)
H, s, out = cyl.zavg([KE2d.transpose()], M.r, M.nr, 1, save=False)

KE_s = out[0]

os.chdir('..')

X = np.vstack((s, KE_s))

np.savetxt('KE_zavg_s.dat', np.transpose(X))

del M
del KE
