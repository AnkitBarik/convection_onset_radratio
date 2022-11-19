import numpy as np
import scipy.sparse as ss
import sys
import numpy.polynomial.chebyshev as ch
from glob import glob
import os

pwd = os.getcwd()

sys.path.insert(0, pwd)
sys.path.insert(0, pwd + '/bin')

import utils as ut
import parameters as par


solnum = 0
nR = par.N + 2# number of radial points

lmax = par.lmax
m    = par.m
symm = par.symm
N    = par.N
Ek   = par.Ek_gap
chi  = par.ricb
ricb = chi/(1-chi)
rcmb = 1/(1-chi)
n    = ut.n

R1 = ricb
R2 = rcmb

# xk are the colocation points, from -1 to 1
i = np.arange(0, nR-2)
xk = np.r_[1, np.cos( (i+0.5)*np.pi/nR), -1]
# rk are the radial colocation points, from Ra to Rb
r = 0.5*(rcmb-ricb)*( xk + 1 ) + ricb
r2 = r**2

Inv_r = ss.diags(1/r, 0)


# matrix with Chebishev polynomials at every x point for all degrees:
chx = ch.chebvander(xk, par.N-1) # this matrix has nR rows and N-1 cols


reu = np.loadtxt('real_flow.field', usecols=solnum)
imu = np.loadtxt('imag_flow.field', usecols=solnum)


Plj0 = reu[:n]    + 1j*imu[:n] 		#  N elements on each l block
Tlj0 = reu[n:n+n] + 1j*imu[n:n+n] 	#  N elements on each l block
Plj   = np.reshape(Plj0, (-1, N))
Tlj   = np.reshape(Tlj0, (-1, N))


# l-indices for u and b:
lmm = 2*np.shape(Plj)[0] -1 # this should be =lmax-m
s = int(symm*0.5+0.5) # s=0 if u is antisymm, s=1 if u is symm
if m>0:
	lup = np.arange( m+1-s, m+1-s +lmm, 2) # u pol
	lut = np.arange( m+s, m+s   +lmm, 2) # u tor
elif m==0:
	lup = np.arange( 1+s, 1+s +lmm, 2) # u pol
	lut = np.arange( 2-s, 2-s +lmm, 2) # u tor


Lup = ss.diags(lup*(lup+1), 0)
Lut = ss.diags(lut*(lut+1), 0)
Inv_Lup = ss.diags(1/(lup*(lup+1)), 0)
Inv_Lut = ss.diags(1/(lut*(lut+1)), 0)


d1Plj = np.zeros(np.shape(Plj), dtype=complex)
d2Plj = np.zeros(np.shape(Plj), dtype=complex)
d1Tlj = np.zeros(np.shape(Tlj), dtype=complex)

P0 = np.zeros((int((lmm+1)/2), nR), dtype=complex)
P1 = np.zeros((int((lmm+1)/2), nR), dtype=complex)
P2 = np.zeros((int((lmm+1)/2), nR), dtype=complex)
T0 = np.zeros((int((lmm+1)/2), nR), dtype=complex)
T1 = np.zeros((int((lmm+1)/2), nR), dtype=complex)

Q0 = np.zeros((int((lmm+1)/2), nR), dtype=complex)
S0 = np.zeros((int((lmm+1)/2), nR), dtype=complex)


np.matmul( Plj, chx.T, P0 )
np.matmul( Tlj, chx.T, T0 )

for k in range(np.size(lup)):
	d1Plj[k,:] = ut.Dcheb(   Plj[k,:], ricb, rcmb)
	d2Plj[k,:] = ut.Dcheb( d1Plj[k,:], ricb, rcmb)
np.matmul(d1Plj, chx.T, P1)
np.matmul(d2Plj, chx.T, P2)

for k in range(np.size(lup)):
	d1Tlj[k,:] = ut.Dcheb(   Tlj[k,:], ricb, rcmb)
np.matmul(d1Tlj, chx.T, T1)


Q0 = Lup * P0 * Inv_r
Q1 = ( Lup * P1 -   Q0 ) * Inv_r

S0 = P1 + P0 * Inv_r
S1 = P2 + Inv_Lup * Q1


dint_pol = np.zeros((int((lmm+1)/2), nR), dtype=complex)
dint_tor = np.zeros((int((lmm+1)/2), nR), dtype=complex)
#dkin_pol = np.zeros((int((lmm+1)/2), nR),dtype=complex)
#dkin_tor = np.zeros((int((lmm+1)/2), nR),dtype=complex)


## Poloidal components
for k, l in enumerate(lup):

	L = l*(l+1)

	q0 = Q0[k,:]
	q1 = Q1[k,:]
	s0 = S0[k,:]
	s1 = S1[k,:]

	f0 = 4*np.pi/(2*l+1)

	# Dint: Internal energy dissipation
	f1 = L*np.absolute(q0 + r*s1 - s0)**2
	f2 = 3*np.absolute(r*q1)**2
	f3 = L*(l-1)*(l+2)*np.absolute(s0)**2
	dint_pol[k,:] = f0*( f1+f2+f3 )

	# # Dkin: Kinetic energy dissipation
	# f1 = L * r2 * np.conj(s0) * s2
	# f2 = 2 * rk * L * np.conj(s0) * s1
	# f3 = -(L**2)*( np.conj(s0)*s0 ) - (l**2+l+2) * ( np.conj(q0)*q0 )
	# f4 = 2 * rk * np.conj(q0)*q1 + r2 * np.conj(q0) * q2
	# f5 = 2 * L *( np.conj(q0)*s0 + q0*np.conj(s0) )
	# dkin_pol[k,:] = f0*( f1+f2+f3+f4+f5 )


## Toroidal components
for k, l in enumerate(lut):

	L = l*(l+1)

	t0 = T0[k,:]
	t1 = T1[k,:]

	f0 = 4*np.pi/(2*l+1)

	# Dint: Internal energy dissipation
	f1 = L*np.absolute( r*t1-t0 )**2
	f2 = L*(l-1)*(l+2)*np.absolute( t0 )**2
	dint_tor[k,:] = f0*( f1+f2 )


dpol = np.real(sum(dint_pol, 0)*par.Ek)
dtor = np.real(sum(dint_tor, 0)*par.Ek)
dint = np.real(sum(dint_pol+dint_tor, 0)*par.Ek)


var = np.vstack([r, dpol, dtor, dint])

np.savetxt('dissipation_profile.dat', var.transpose())
