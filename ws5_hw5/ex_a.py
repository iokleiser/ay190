#!/usr/bin/env python

import os
import array
import math
from numpy import *
import scipy
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.figure as fig

# finite differencing
# (i) forward differencing

def find_root_bisection(x_up,x_dn,func,tol):
    frac_err = (x_up - x_dn)/x_dn
    x_mid = (x_up + x_dn)/2
    iterations = 0
    while frac_err >= tol:
        #print "iteration: ", iterations
        #print "upper x, lower x: ",x_up,", ",x_dn
        #print "fractional error: ",frac_err
        f_up = func(x_up)
        f_dn = func(x_dn)
        if f_up * f_dn <= 0:
            x_mid = (x_up + x_dn)/2.0
            f_mid = func(x_mid)
        else:
            print "problem!! x values do not bracket a root"
            quit()
        if f_up * f_mid < 0:
            x_dn = x_mid
        elif f_dn * f_mid < 0:
            x_up = x_mid
        frac_err = (x_up - x_dn)/x_dn
        iterations += 1
    return x_mid,iterations
            




# part 1
# values we're given:
days_to_seconds = 24 * 60 * 60
period = 365.25635 * days_to_seconds # earth period in seconds
omega = 2*math.pi/period
a = 1.496e6
ep = 0.0167
b = a*sqrt(1-ep**2)
time = array([91, 182, 273]) *  days_to_seconds
tol = 1e-10 # fractional error allowed

E_up = 2
E_dn = 1
t = time[0]
eccentric_anomaly = lambda x: x - omega * t - ep * sin(x)

E0,it0 = find_root_bisection(E_up,E_dn,eccentric_anomaly,tol)

E_up = 4
E_dn = 2
t = time[1]
eccentric_anomaly = lambda x: x - omega * t - ep * sin(x)

E1,it1 = find_root_bisection(E_up,E_dn,eccentric_anomaly,tol)

E_up = 6.6
E_dn = 4
t = time[2]
eccentric_anomaly = lambda x: x - omega * t - ep * sin(x)

E2,it2 = find_root_bisection(E_up,E_dn,eccentric_anomaly,tol)

print "eccentricity: ",ep
print "t = ",time[0],", E = ",E0,", iterations = ",it0
print "x = ",a*cos(E0),", y = ",b*sin(E0)
print "t = ",time[1],", E = ",E1,", iterations = ",it1
print "x = ",a*cos(E1),", y = ",b*sin(E1)
print "t = ",time[2],", E = ",E2,", iterations = ",it2
print "x = ",a*cos(E2),", y = ",b*sin(E2)

ep = 0.99999

E_up = 4
E_dn = 1
t = time[0]
eccentric_anomaly = lambda x: x - omega * t - ep * sin(x)

E0,it0 = find_root_bisection(E_up,E_dn,eccentric_anomaly,tol)

E_up = 4
E_dn = 1
t = time[1]
eccentric_anomaly = lambda x: x - omega * t - ep * sin(x)

E1,it1 = find_root_bisection(E_up,E_dn,eccentric_anomaly,tol)

E_up = 4
E_dn = 1
t = time[2]
eccentric_anomaly = lambda x: x - omega * t - ep * sin(x)

E2,it2 = find_root_bisection(E_up,E_dn,eccentric_anomaly,tol)
print "eccentricity: ",ep
print "t = ",time[0],", E = ",E0,", iterations = ",it0
print "x = ",a*cos(E0),", y = ",b*sin(E0)
print "t = ",time[1],", E = ",E1,", iterations = ",it1
print "x = ",a*cos(E1),", y = ",b*sin(E1)
print "t = ",time[2],", E = ",E2,", iterations = ",it2
print "x = ",a*cos(E2),", y = ",b*sin(E2)
