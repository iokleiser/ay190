#!/usr/bin/env python

import sys,math
from pylab import *
from numpy import *
from math import *

def apply_bcs(x,y):
    # apply boundary conditions
    # you need to fill in code
    y[0] = 0
    y[-1] = 0
    return y


def analytic(xdata,x0,sigma):
    func = lambda x: exp(-(x - x0)**2/(2*sigma**2))
    y = zeros(len(x))
    for i in range(len(x)):
        y[i] = func(x[i])
    return y

def upwind(x,y,dx,v,dt):
    length = len(y)
    ynew = array(zeros(length))
    ynew[1:] = y[1:] - v*dt/dx * (y[1:] - y[:(length-1)])
    return ynew

def lax_friedrich(x,y,dx,v,dt):
    length = len(y)
    ynew = array(zeros(length))
    ynew[2:] = 0.5*(y[2:] - y[:(length-2)]) - v*dt/dx * (y[2:] - y[:(length-2)])
    return ynew

def ftcs(x,y,dx,v,dt):
    length = len(y)
    ynew = array(zeros(length))
    ynew[2:] = y[1:(length-1)] - v*dt/dx * (y[2:] - y[:(length-2)])
    return ynew

def leapfrog(x,y,dx,v,dt):

# parameters
dx = 0.1
v = 0.1
# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x = arange(0,100,dx)
n = len(x)
y = zeros(n)
cfl = 1.0
dt = 0.5
t = 0.0

# for initial data
sigma = sqrt(15.0)
x0 = 30.0

#set up initial conditions
y = analytic(x,x0,sigma)
print y
print len(x)
print len(y)
# evolve (and show evolution)
ion()
figure()
plot(x,y,'x-') # numerical data
plot(x,analytic(x,x0,sigma),'r-') # analytic data
show()

yold2 = y
yold = y
ntmax = 500

for it in range(ntmax):
    t += dt
    # save previous and previous previous data
    yold2 = yold
    yold = y

    # get new data; ideally just call a function
    y = upwind(x,yold,dx,v,dt)

    # after update, apply boundary conditions
    apply_bcs(x,y) 
    #[FILL IN CODE]

    # get analytic result for time t
    x0 += v * dt
    yana = analytic(x,x0,sigma)
    # compute error estimage
    err = mean(abs(y - yana))
    print "it = ",it,err
    clf()
    # plot numerical result
    plot(x,y,'k')
    # plot analytic results
    plot(x,yana,'r-')
    draw()


show()


