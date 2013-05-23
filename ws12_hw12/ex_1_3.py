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


pc = genfromtxt('./pc.dat')
minmod = genfromtxt('./minmod.dat')
mc = genfromtxt('./mc.dat')
x = pc[:,0]
x1 = x[x < 2.4180111E-01]
x3a = x[x > 4.6319762E-01]
x3 = x3a[x3a < 6.6604739E-01]
#x3 = x[x > 4.6319762E-01 and x < 6.6604739E-01]
x4a = x[x > 6.6604739E-01]
x4 = x4a[x4a < 9.1987044E-01]
#x4 = x[x > 6.6604739E-01 and x < 9.1987044E-01]
x5 = x[x > 9.1987044E-01]

y1 = 1.0 + zeros(len(x1))
y3 = 4.8490930E-01 + zeros(len(x3))
y4 = 1.6541857E-01 + zeros(len(x4))
y5 = 1.0000000E-01 + zeros(len(x5))

x = concatenate((x1,x3,x4,x5))
y = concatenate((y1,y3,y4,y5))

plt.clf()
plt.plot(pc[:,0],pc[:,1],linewidth=2)
plt.plot(minmod[:,0],minmod[:,1],'r',linewidth=2)
plt.plot(mc[:,0],mc[:,1],'k',linewidth=2)
plt.plot(x,y,'k--',linewidth=2)
plt.xlabel('x')
plt.ylabel('$\\rho$')
leg = plt.legend(['pc','minmod','mc','analytical'],loc=0)
leg.draw_frame(False)
plt.savefig('profiles.pdf',format='pdf')

