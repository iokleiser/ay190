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


def analytic(x,x0,sigma):
    y = lambda x: exp(-(x - x0)**2/(2*sigma**2))

def upwind(x,y,dx,v,dt):
    ynew = array(zeros(len(y)))
    for j in arange(1,len(y))
        ynew[j] = y[j] - v*dt/dx * (y[j] - y[j-1])
    return ynew

def downwind(x,y,dx,v,dt):
    ynew = array(zeros(len(y)))
    for j in arange(0,len(y)-1)
        ynew[j] = y[j] - v*dt/dx * (y[j+1] - y[j])
    return ynew

# parameters
dx = 0.1
v = 0.1
# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x = arange(0,100,dx)
n = len(x)
y = zeros(n)
cfl = 1.0
dt = 0.1
t = 0.0

# for initial data
sigma = sqrt(15.0)
x0 = 30.0

#set up initial conditions
y = analytic(x,x0,sigma)

# evolve (and show evolution)
ion()
figure()
plot(x,y,'x-') # numerical data
plot(x,analytic(x,x0,sigma),'r-') # analytic data
show()

yold2 = y
yold = y
ntmax = 2000
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
    err = 0
    # err = ???
    #[FILL IN CODE]
    print "it = ",it,err
    clf()
    # plot numerical result
    plot(x,y,'k')
    # plot analytic results
    plot(x,yana,'r-')
    draw()


show()


