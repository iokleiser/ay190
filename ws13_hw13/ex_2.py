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
import numpy.random as random

def randwalk_1D(lam,xr,x0,s):
    #random.seed(s)
    stp = 0
    x = x0
    xvec = []
    count = 0
    counts = []
    while stp == 0:
        r = random.rand() - 0.5
        if r >= 0:
            step = lam
        elif r < 0:
            step = -lam
        x += step
        xvec.append(x)
        if x <= xr[0] or x >= xr[1]:
            stp = 1
        count += 1
        counts.append(count)
    return xvec,counts

def randwalk_2D(lam,xr,yr,r0,s):
    directions = [[0,1],[1,0],[-1,0],[0,-1],[1.0/sqrt(2),1.0/sqrt(2)],[1.0/sqrt(2),-1.0/sqrt(2)],[-1.0/sqrt(2),-1.0/sqrt(2)],[-1.0/sqrt(2),1.0/sqrt(2)]]
    #random.seed(s)
    stp = 0
    x = r0[0]
    y = r0[1]
    xvec = []
    yvec = []
    count = 0
    counts = []
    while stp == 0:
        ind = int(floor(random.rand() * 8))
        d = directions[ind]
        x += lam * d[0]
        y += lam * d[1]
        xvec.append(x)
        yvec.append(y)
        if x <= xr[0] or x >= xr[1] or y <= yr[0] or y >= yr[1]:
            stp = 1
        count += 1
        counts.append(count)
    return xvec,yvec,counts


xr = [-1.0,1.0]
lam = 0.01
x0 = 0
s = 100
random.seed(s)
xvec,counts = randwalk_1D(lam,xr,x0,s)
plt.xlabel('counts')
plt.ylabel('x')
plt.plot(counts,xvec)
plt.savefig('randwalk_1D.pdf',format='pdf')

yr = [-1.0,1.0]
r0 = [0.0,0.0]
random.seed(1)
x1,y1,counts1 = randwalk_2D(lam,xr,yr,r0,s)
random.seed(10)
x2,y2,counts2 = randwalk_2D(lam,xr,yr,r0,s)
random.seed(100)
x3,y3,counts3 = randwalk_2D(lam,xr,yr,r0,s)
plt.clf()
plt.xlabel('x')
plt.ylabel('y')
plt.plot(x1,y1)
plt.plot(x2,y2)
plt.plot(x3,y3)
plt.savefig('randwalk_2D.pdf',format='pdf')
