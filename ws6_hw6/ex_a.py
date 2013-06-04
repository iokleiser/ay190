#!/usr/bin/env python

import os
import array
import math
from numpy import *
import scipy
from scipy.linalg import solve
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.figure as fig

def gauss_elim(m,b):
    ncol = m.shape[1]
    nrow = m.shape[0]
    b = reshape(b,(len(b),1))
    totmx = hstack((m,b))
    triangular = 0
    while triangular == 0:
        for i in range(nrow):
            a = totmx[i,i]
            for j in arange(i + 1,nrow):
                irow = totmx[i,:]
                jrow = totmx[j,:]
                newrow = -1*totmx[j,i]/a * irow + jrow
                totmx[j,:] = newrow
            triangular = check_if_triang(totmx)
    
    return totmx

def check_if_triang(m):
    nrow = m.shape[0]
    triangular = []
    for i in range(nrow):
        row = m[i,:]
        triangular.append(all(row[0:i+1] == 0))
    return all(triangular)


LSE1_m = genfromtxt('LSE1_m.dat')
LSE1_bvec = genfromtxt('LSE1_bvec.dat')

LSE1_g = gauss_elim(LSE1_m,LSE1_bvec)
print LSE1_g

x = solve(LSE1_m,LSE1_bvec)


