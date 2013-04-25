#!/opt/local/bin/python

import sys
from scipy import *
from pylab import *
from math import *

# global constants
ggrav = 6.67e-8
clite = 3.0e10
msun = 1.99e33


# EOS
# neutron stars:
# polyG = 2.0
# polyK = 100.0 * 5.55e38/6.1755e17**polyG

# EOS for 
# white dwarfs:
polyG = 4.0/3.0
polyK = 1.244e15*0.5**polyG

# central values
rhoc = 1.0e10

# minimum pressure
rhomin = 1.0e-10*rhoc
min_press = polyK * rhomin**polyG


# grid
rmax = 1.0e9


def set_grid(rmax,nzones):
    # set up the grid and return the
    # radius array and dr

    rad = array(arange(0.0,rmax,nzones))
    dr = rad[1] - rad[0]

    return (rad,dr)

def tov_RHS(r,data):

    rhs = zeros(2)

    mass = data[1]
    press = max(data[0],min_press)

    rho = (press/polyK)**(1/polyG)

    rhs[0] = -1*ggrav*mass*rho/r**2
    rhs[1] = 4*pi*rho*r**2

    return rhs

def tov_RK2(old_data,r,dr):

    [ FILL IN RK2 ]

    return new_data
    
def tov_RK3(old_data,r,dr):

    [ FILL IN RK3 ]

    return new_data


def tov_RK4(old_data,r,dr):

    [ FILL IN RK4 ]

    return new_data


def tov_integrate(rmax,nzones):

    # set up grid
    (rad,dr) = set_grid(rmax,nzones)

    # initialize some variables
    tovdata = zeros((nzones,2))
    # 0 -- press
    # 1 -- mbary

    tovout = zeros((nzones,4))
    # 0 -- rho
    # 1 -- press
    # 2 -- eps
    # 3 -- mass
     
    # central values
    tovdata[0,0] = polyK * rhoc**polyG
    tovdata[0,1] = 0.0


    # you will need to track the surface (where press <= press_min)
    isurf = 0
    found_surface = 0
    for i in range(nzones-1):

        # integrate one step using RK2 (can change this to
        # RK3 or RK4)
        tovdata[i+1,:] = tov_RK2(tovdata[i,:],rad[i],dr)

        # check if press below 0
        if(tovdata[i+1,0] <= min_press and found_surface == 0):
            isurf = i
            found_surface = 1

        # press and mass
        tovout[i+1,1] = tovdata[i+1,0]

        if (i+1 > isurf and isurf > 0):
            tovout[i+1,3] = tovdata[isurf,1]
        else:
            tovout[i+1,3] = tovdata[i+1,1]
        # compute density
        tovout[i+1,0] = (tovout[i+1,1]/polyK)**(1/polyG)
        # compute eps
        tovout[i+1,2] = tovout[i+1,1]/((polyG - 1)*tovout[i+1,0])

    return (tovout,isurf,dr)

# for convergence: 
# number of points
na = array([10.0,20.0,50.0,100.0,200.0,500.0,1.0e3,2.0e3,5.0e3,1.0e4] ])
# to store masses
masses = zeros(len(na))
# to store the drs
drs = zeros(len(na))

for i in range(len(na)):
    (tov_star,isurf,dr) = tov_integrate(rmax,na[i])
    [ FILL IN CODE FOR CONVERGENCE CALCULATION ]
