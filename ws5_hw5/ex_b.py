#!/usr/bin/env python

import os
import array
import math
from numpy import *
import scipy
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.figure as fig


def find_root_newton(init_guess,func,epsilon,tol):
    x = init_guess
    iterations = 0
    dx = tol * 1.0e5
    while abs(dx) > tol:
        xf = x*(1 + epsilon)
        xb = x*(1 - epsilon)
        fprime = (func(xf) - func(xb))/(xf - xb)
        dx = -1*func(x)/fprime
        print x
        x += dx
        iterations += 1
        if iterations > 300:
            iterations = 0
            print "did not converge!"
            conv = 0
            return x,iterations,conv
    conv = 1
    return x,iterations,conv

def multi_newton(init_grid,func,epsilon,tol):
    conv = 0.0
    i = 0
    while conv == 0:
        init_x = init_grid[i]
        xr,iterations,conv = find_root_newton(init_x,func,epsilon,tol)
        if conv == 0:
            i += 1
    if conv == 1:
        newfunc = lambda x: func/(x - xr)

func = lambda x: 3*x**5 + 5*x**4 - x**3

x,iterations = find_root_newton(1,func,0.01,1.0e-10)
print x
print iterations
