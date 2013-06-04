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

def calc_pi(N):
    random.seed(1)
    Nin = 0
    for i in range(N):
        x = random.rand()
        y = random.rand()
        r = sqrt(x**2 + y**2)
        if r <= 1:
            Nin += 1
    pival = 4.0*Nin/N
    return pival



Nlist = [100,200,500,1000,2000,5000,1e4,2e4,5e4]
pival = []
for N in Nlist:
    newpi = calc_pi(int(N))
    pival.append(newpi)
Nlist = array(Nlist)
pival = array(pival)
plt.clf()
plt.plot(Nlist,abs(pival - pi),'r',linewidth=2)
plt.savefig('pi_mc.pdf',format='pdf')




