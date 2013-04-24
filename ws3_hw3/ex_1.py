#!/usr/bin/env python

import os
import array
import math
from numpy import *
import scipy
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.figure as fig

# part a

def trap_integ(fn,xr,h):
   x = array(arange(xr[0],xr[1],h))
   n = len(x)
   Qi = []
   for i in range(n-1):
       a = x[i]
       b = x[i+1]
       fa = fn(a)
       fb = fn(b)
       Qi.append((1.0/2.0)*(b - a)*(fb + fa))
   Q = sum(Qi)
   return Q



def simpson_integ(fn,xr,h):
   x = array(arange(xr[0],xr[1],h))
   n = len(x)
   Qi = []
   for i in range(n-1):
       a = x[i]
       b = x[i+1]
       fa = fn(a)
       fb = fn(b)
       fm = fn((a+b)/2)
       Qi.append((b - a)/6 * (fa + 4*fm + fb))
   Q = sum(Qi)
   return Q

fn = lambda x: sin(x)

hvec = arange(0.0001,0.1,0.0001)
xr = [0,math.pi]


Qtrap = []
Qsimp = []
for h in hvec:
   Qtrap.append(trap_integ(fn,xr,h))
   Qsimp.append(simpson_integ(fn,xr,h))


plt.clf()
plt.plot(hvec,Qtrap)
plt.plot(hvec,Qsimp)
plt.axis([0.0,0.10,1.993,2.001])
plt.xlabel("stepsize h")
plt.ylabel("integral of sin(x) on [0,$\pi$]")
lg = plt.legend(["Trapezoidal","Simpson"],loc=0)
lg.draw_frame(False)
plt.savefig("trapsimp_sinx.pdf",format="pdf")


fn = lambda x: x*sin(x)

Qtrap = []
Qsimp = []
for h in hvec:
   Qtrap.append(trap_integ(fn,xr,h))
   Qsimp.append(simpson_integ(fn,xr,h))

plt.clf()
plt.plot(hvec,Qtrap)
plt.plot(hvec,Qsimp)
#plt.axis([0.0,0.10,1.993,2.001])
plt.xlabel("stepsize h")
plt.ylabel("integral of xsin(x) on [0,$\pi$]")
lg = plt.legend(["Trapezoidal","Simpson"],loc=0)
lg.draw_frame(False)
plt.savefig("trapsimp_xsinx.pdf",format="pdf")
