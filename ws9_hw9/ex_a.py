#!/usr/bin/env python

import os
import array
from math import *
from numpy import *
from scipy import *
from scipy.interpolate import *
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.figure as fig


arr = genfromtxt('./presupernova.dat')

mass = arr[:,1]
radius = arr[:,2]
rho = arr[:,4]

#1
plt.clf()
plt.semilogy(radius,rho,'r',linewidth=2)
plt.xlabel('radius (cm)')
plt.ylabel('density')
plt.xlim([1e9,1e13])

plt.savefig('rho.pdf',format='pdf')


def make_grid(radius,rho,dx,xmax):
    
    x = arange(dx,xmax,dx)

    rhoSpline = InterpolatedUnivariateSpline(radius,rho)
    rhoNew = rhoSpline(x)
    return x,rhoNew

def get_Phi(phiOut,rhoNew,x,dx):
    z = zeros(len(x)) #dPhi/dx
    Phi = zeros(len(x))
    for i in range(len(x)-1):
        z[i+1] = z[i] + dx * (4*pi*G*rhoNew[i] - 2*z[i]/x[i])
        Phi[i+1] = Phi[i] + dx*z[i]
    
    Phi += phiOut - Phi[-1]
    return Phi

#2
dx = 1.0e6
xmax = 1.0e9

x,rhoNew = make_grid(radius,rho,dx,xmax)


#3
M = sum(4*pi*dx*rhoNew*x**2)
G = 6.67e-8
phiOut = -1*G*M/xmax
Phi = get_Phi(phiOut,rhoNew,x,dx)

plt.clf()
plt.plot(x,Phi,linewidth=2)
plt.xlabel('radius')
plt.ylabel('$\\Phi$')
plt.savefig('phi.pdf',format='pdf')

vol = 4.0*pi/3.0*xmax**3
rhoNew = ones(len(x))*M/vol
Phi = get_Phi(phiOut,rhoNew,x,dx)
smallphi = 2.0/3.0 * pi * G * rhoNew * (x**2 - 3*(x[-1])**2)

plt.clf()
plt.plot(x,Phi,linewidth=2)
plt.plot(x,smallphi,'r',linewidth=2)
plt.xlabel('radius')
plt.ylabel('$\\Phi$')
leg = plt.legend(['Euler','analytical'],loc=0)
leg.draw_frame(False)
plt.savefig('phi2.pdf',format='pdf')

dxs = [1.0e5, 2.0e5, 5.0e5, 1.0e6, 2.0e6, 5.0e6, 1.0e7, 2.0e7, 5.0e7]
rel_errs = []
for dx in dxs:
    Phi = get_Phi(phiOut,rhoNew,x,dx)
    diff = Phi[-1] - smallphi[-1]
    err = mean(abs(Phi-smallphi-diff))
    print err
    rel_errs.append(err)
rel_errs = array(rel_errs)
plt.clf()
plt.xlabel('radius')
plt.ylabel('dx')
plt.loglog(dxs,rel_errs,'kx',linewidth=2)
plt.savefig('rel_err.pdf',format='pdf')
