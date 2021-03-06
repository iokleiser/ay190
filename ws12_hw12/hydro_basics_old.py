#!/usr/bin/env python
import sys,math
from pylab import *

# basic parameters
gamma = 5.0/3.0
cfl = 0.5
dt = 1.0e-5
dtp = dt
reconstruction_type = 'mc' # minmod, mc, pc
# use nzones
nzones = 200
# run until time 0.2
tend = 0.2

##################################### class definition
class mydata:
    def __init__(self,nzones):
        self.x     = zeros(nzones) # cell centers
        self.xi    = zeros(nzones) # cell LEFT interfaces
        self.rho   = zeros(nzones)
        self.rhop   = zeros(nzones)
        self.rhom   = zeros(nzones)
        self.vel    = zeros(nzones)
        self.velp   = zeros(nzones)
        self.velm   = zeros(nzones)
        self.eps    = zeros(nzones)
        self.epsp   = zeros(nzones)
        self.epsm   = zeros(nzones)
        self.press  = zeros(nzones)
        self.pressp = zeros(nzones)
        self.pressm = zeros(nzones)
        self.q     = zeros((3,nzones)) # conserved quantities
        self.qp    = zeros((3,nzones))  
        self.qm    = zeros((3,nzones))  
        self.n     = nzones
        self.g     = 3

    def setup_grid(self,xmin,xmax):
        dx = (xmax - xmin) / (self.n - self.g*2 - 1)
        xmin = xmin - self.g*dx
        xmax = xmax + self.g*dx
        for i in range(self.n):
            # cell centers
            self.x[i] = xmin + (i)*dx
        # cell LEFT interfaces
        for i in range(self.n):
            self.xi[i] = self.x[i] - 0.5*dx

    
    def setup_ID(self):
        # Shocktube initial data
        rchange = (self.x[self.n-self.g-1] - self.x[self.g]) / 2.0
        print rchange
        rho1 = 1.0
        rho2 = 0.1
        press1 = 1.0
        press2 = 0.125
        for i in range(self.n):
            if self.x[i] < rchange:
                self.rho[i] = rho1
                self.press[i] = press1
                self.eps[i] = press1 / (rho1) / (gamma - 1.0)
                self.vel[i] = 0.0
            else:
                self.rho[i] = rho2
                self.press[i] = press2
                self.eps[i] = press2 / (rho2) / (gamma - 1.0)
                self.vel[i] = 0.0


##################################### some basic functions
def prim2con(rho,vel,eps):
    q = zeros((3,len(rho)))
    q[0,:] = rho[:]
    q[1,:] = rho[:] * vel[:]
    q[2,:] = rho[:] * eps[:] + \
                  0.5 * rho[:] * vel[:]**2

    return q

def con2prim(q):
    rho = q[0,:]
    vel = q[1,:] / rho
    eps = q[2,:] / rho - 0.5*vel**2
    press = eos_press(rho,eps,gamma)

    return (rho,eps,press,vel)

############# boundary conditions
def apply_bcs(hyd):
    hyd.rho[0:hyd.g-1] = hyd.rho[hyd.g]
    hyd.vel[0:hyd.g-1] = hyd.vel[hyd.g]
    hyd.eps[0:hyd.g-1] = hyd.eps[hyd.g]
    hyd.press[0:hyd.g-1] = hyd.press[hyd.g]

    hyd.rho[hyd.n-hyd.g:hyd.n-1] = hyd.rho[hyd.n-hyd.g-1]
    hyd.vel[hyd.n-hyd.g:hyd.n-1] = hyd.vel[hyd.n-hyd.g-1]
    hyd.eps[hyd.n-hyd.g:hyd.n-1] = hyd.eps[hyd.n-hyd.g-1]
    hyd.press[hyd.n-hyd.g:hyd.n-1] = hyd.press[hyd.n-hyd.g-1]

    return hyd

############# minmod function 
def minmod(a,b):
    if(a*b < 0):
        mm = 0.0
    elif(abs(a)<abs(b)):
        mm=a
    else:
        mm=b
    return mm

############# minmod function 
def tvd_minmod_reconstruct(n,g,f,x,xi):
    fp = zeros(n)
    fm = zeros(n)
    for i in range(g-1,n-g+1):
        df_up = (f[i] - f[i-1])/(x[i] - x[i-1])
        df_down = (f[i+1] - f[i])/(x[i+1] - x[i])
        Delta = minmod(df_up,df_down)
        fm[i] = f[i] - Delta * (x[i] - xi[i])
        fp[i] = f[i] + Delta * (xi[i+1] - x[i])

    return (fp,fm)

############# signum functions
def signum(x,y):
    if(y >= 0):
        return abs(x)
    else:
        return -abs(x)

############# mc reconstruction
def tvd_mc_reconstruct(n,g,f,x,xi):
    fp = zeros(n)
    fm = zeros(n)
    for i in range(g-1,n-g+1):
        df_up = (f[i] - f[i-1])/(x[i] - x[i-1])
        df_down = (f[i+1] - f[i])/(x[i+1] - x[i])
        if (df_up * df_down < 0):
            Delta = 0
        else:
            Delta = signum(min(2*abs(df_up),2*abs(df_down),0.5*(abs(df_up) + abs(df_down))), (df_up + df_down))
        fm[i] = f[i] - Delta * (x[i] - xi[i])
        fp[i] = f[i] + Delta * (xi[i+1] - x[i])

    return (fp,fm)


############# reconstruction top level function
def reconstruct(hyd,type):
    if(type=='pc'):
        # piecewise constant reconstruction 
        for i in range(hyd.g-1,hyd.n-hyd.g+1):
            hyd.rhop[i] = hyd.rho[i]
            hyd.rhom[i] = hyd.rho[i]
            hyd.epsp[i] = hyd.eps[i]
            hyd.epsm[i] = hyd.eps[i]
            hyd.velp[i] = hyd.vel[i]
            hyd.velm[i] = hyd.vel[i]


    elif(type=='minmod'):
        (hyd.rhop,hyd.rhom) = tvd_minmod_reconstruct(hyd.n,hyd.g,hyd.rho,hyd.x,hyd.xi)
        (hyd.epsp,hyd.epsm) = tvd_minmod_reconstruct(hyd.n,hyd.g,hyd.eps,hyd.x,hyd.xi)
        (hyd.velp,hyd.velm) = tvd_minmod_reconstruct(hyd.n,hyd.g,hyd.vel,hyd.x,hyd.xi)

    elif(type=='mc'):
        (hyd.rhop,hyd.rhom) = tvd_mc_reconstruct(hyd.n,hyd.g,hyd.rho,hyd.x,hyd.xi)
        (hyd.epsp,hyd.epsm) = tvd_mc_reconstruct(hyd.n,hyd.g,hyd.eps,hyd.x,hyd.xi)
        (hyd.velp,hyd.velm) = tvd_mc_reconstruct(hyd.n,hyd.g,hyd.vel,hyd.x,hyd.xi)
                
    else:
        print "reconstruction type not known; abort!"
        sys.exit()


    hyd.pressp = eos_press(hyd.rhop,hyd.epsp,gamma)
    hyd.pressm = eos_press(hyd.rhom,hyd.epsm,gamma)

    hyd.qp = prim2con(hyd.rhop,hyd.velp,hyd.epsp)
    hyd.qm = prim2con(hyd.rhom,hyd.velm,hyd.epsm)

    return hyd

############# equation of state
def eos_press(rho,eps,gamma):
    press = (gamma - 1.0) * rho * eps
    return press

def eos_cs2(rho,eps,gamma):
    prs = (gamma - 1.0) * rho *eps
    dpde = (gamma - 1.0) * rho
    dpdrho = (gamma - 1.0) * eps
    cs2 = dpdrho + dpde * prs/(rho+1.0e-30)**2
    return abs(cs2)

############# time step calculation
def calc_dt(hyd,dtp):
    cs = sqrt(eos_cs2(hyd.rho,hyd.eps,gamma))
    dtnew = 1.0
    for i in range(hyd.g,hyd.n-hyd.g):
        dtnew = min(dtnew, (hyd.x[i+1]-hyd.x[i]) / \
                    max(abs(hyd.vel[i]+cs[i]), \
                        abs(hyd.vel[i]-cs[i])))

    dtnew = min(cfl*dtnew,1.05*dtp)
    return dtnew

############# HLLE solver
def hlle(hyd):
    fluxdiff = zeros((3,hyd.n))
    # compute eigenvalues
    evl  = zeros((3,hyd.n))
    evr  = zeros((3,hyd.n))
    smin = zeros(hyd.n)
    smax = zeros(hyd.n)
    csp  = sqrt(eos_cs2(hyd.rhop,hyd.epsp,gamma))
    csm  = sqrt(eos_cs2(hyd.rhom,hyd.epsm,gamma))
    for i in range(1,hyd.n-2):
        evl[0,i] = hyd.velp[i]
        evl[1,i] = hyd.velp[i] - csp[i]
        evl[2,i] = hyd.velp[i] + csp[i]
        evr[0,i] = hyd.velm[i+1]
        evr[1,i] = hyd.velm[i+1] - csm[i+1]
        evr[2,i] = hyd.velm[i+1] + csm[i+1]
        smin[i] = min(evl[0,i],evl[1,i],evl[2,i],evr[0,i],evr[1,i],evr[2,i],0.0)
        smax[i] = max(evl[0,i],evl[1,i],evl[2,i],evr[0,i],evr[1,i],evr[2,i],0.0)

    # set up flux left L and right R of the interface
    # at i+1/2
    fluxl = zeros((3,hyd.n))
    fluxr = zeros((3,hyd.n))


    for i in range(1,hyd.n-2):
        fluxl[0,i] = hyd.qp[0,i] * hyd.velp[i]
        fluxl[1,i] = hyd.qp[1,i] * hyd.velp[i] + hyd.pressp[i]
        fluxl[2,i] = (hyd.qp[2,i] + hyd.pressp[i]) * hyd.velp[i]
        fluxr[0,i] = hyd.qm[0,i+1] * hyd.velm[i+1]
        fluxr[1,i] = hyd.qm[1,i+1] * hyd.velm[i+1] + hyd.pressm[i+1]
        fluxr[2,i] = (hyd.qm[2,i+1] + hyd.pressm[i+1]) * hyd.velm[i+1]

    # solve the Riemann problem for the i+1/2 interface
    ds = smax - smin
    flux = zeros((3,hyd.n))
    for i in range(hyd.g-1,hyd.n-hyd.g+1):
        flux[:,i] = (smax[i] * fluxl[:,i] - smin[i] * fluxr[:,i] + smin[i] * smax[i] * (hyd.qm[:,i+1] - hyd.qp[:,i])) / (smax[i] - smin[i])

    # flux differences
    for i in range(hyd.g,hyd.n-hyd.g):
        dx = hyd.xi[i+1] - hyd.xi[i]
        fluxdiff[:,i] = 1.0/dx * (flux[:,i] - flux[:,i-1])

    return fluxdiff

############# RHS calculation
def calc_rhs(hyd):
    # reconstruction and prim2con
    hyd = reconstruct(hyd,reconstruction_type)
    # compute flux differences
    fluxdiff = hlle(hyd)
    # return RHS = - fluxdiff
    return -fluxdiff

########################################################################
# Main program
########################################################################

# initialize
hyd = mydata(nzones)

# set up grid
hyd.setup_grid(0.0,1.0)

# set up initial data
hyd.setup_ID()

# get initial timestep
dt = calc_dt(hyd,dt)

# initial prim2con
hyd.q = prim2con(hyd.rho,hyd.vel,hyd.eps)

t = 0.0
i = 0

# display stuff
ion()
figure()
plot(hyd.x,hyd.rho,"r-")
show()

# main integration loop
while(t < tend):

    if(i % 10 == 0):
        # output
        print "%5d %15.6E %15.6E" % (i,t,dt)
        clf()
        plot(hyd.x,hyd.rho,"r-")
        draw()

    # calculate new timestep
    dt = calc_dt(hyd,dt)

    # save old state
    hydold = hyd
    qold = hyd.q

    # calc rhs
    k1 = calc_rhs(hyd)
    # calculate intermediate step
    hyd.q = qold + 1.0/2.0 * dt * k1
    # con2prim
    (hyd.rho,hyd.eps,hyd.press,hyd.vel) = con2prim(hyd.q)
    # boundaries
    hyd = apply_bcs(hyd)

    #calc rhs
    k2 = calc_rhs(hyd)
    #apply update
    hyd.q = qold + dt * (0.5 * k1 + 0.5 * k2)
    # con2prim
    (hyd.rho,hyd.eps,hyd.press,hyd.vel) = con2prim(hyd.q)
    # apply bcs
    hyd = apply_bcs(hyd)

    # update time
    t = t + dt
    i = i + 1

# display final result
ioff()
clf()
plot(hyd.x,hyd.rho,"r-")
xlabel('x')
ylabel('$\\rho$')
savefig("profile.pdf", format="pdf")

f = open("mc.dat","w")
for i in range(len(hyd.x)):
    f.write("%15.6E %15.6E \n" % (hyd.x[i],hyd.rho[i]))
f.close()
        


