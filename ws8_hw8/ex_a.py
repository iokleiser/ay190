#!/usr/bin/env python

import os
import array
import math
from numpy import *
import scipy
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.figure as fig

def readmod(filename):
    modfile = open(filename, "r")
    linesCounter = 0
    arr = []
    for line in modfile:
        if linesCounter == 2:
            colNames = line.split()
        if linesCounter > 2:
            dataline = map(float,line.split())
            arr.append(dataline)
        linesCounter += 1
    arr = array(arr)
    print colNames
    return arr,colNames

#arr0,col0 = readmod('pp123/myburn_0000_z1.dat')
#arr1,col1 = readmod('pp123/myburn_0001_z1.dat')
arr2,col2 = readmod('pp123/myburn_0002_z1.dat')

plt.clf()
t = arr2[:,col2.index('time')]/(3.15e7)
h1 = arr2[:,col2.index('h1')]
he4 = arr2[:,col2.index('he4')]
plt.semilogy(t,h1,'r',linewidth=2)
plt.semilogy(t,he4,'b',linewidth=2)
plt.xlabel('time (yr)')
plt.ylabel('mass fraction')
leg = plt.legend(['h1','he3'],loc=0)
leg.draw_frame(False)
plt.savefig('T_3e7.pdf',format='pdf')
