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
nxzones = 200
nyzones = 100
# run until time 0.2
tend = 0.2

##################################### class definition
class mydata:
    def __init__(self,nxzones,nyzones):
        self.x     = zeros(nxzones) # cell centers
        self.y     = zeros(nyzones)
        self.xi    = zeros(nxzones) # cell LEFT interfaces
        self.yi    = zeros(nyzones)
        self.rho   = zeros((nxzones,nyzones))
        self.rhopx   = zeros((nxzones,nyzones))
        self.rhopy   = zeros((nxzones,nyzones))
        self.rhomx   = zeros((nxzones,nyzones))
        self.rhomy   = zeros((nxzones,nyzones))
        self.vel    = zeros((nxzones,nyzones))
        self.velpx   = zeros((nxzones,nyzones))
        self.velpy   = zeros((nxzones,nyzones))
        self.velmx   = zeros((nxzones,nyzones))
        self.velmy   = zeros((nxzones,nyzones))
        self.eps    = zeros((nxzones,nyzones))
        self.epspx   = zeros((nxzones,nyzones))
        self.epspy   = zeros((nxzones,nyzones))
        self.epsmx   = zeros((nxzones,nyzones))
        self.epsmy   = zeros((nxzones,nyzones))
        self.press  = zeros((nxzones,nyzones))
        self.presspx = zeros((nxzones,nyzones))
        self.presspy = zeros((nxzones,nyzones))
        self.pressmx = zeros((nxzones,nyzones))
        self.pressmy = zeros((nxzones,nyzones))
        self.q     = zeros((3,nxzones,nyzones)) # conserved quantities
        self.qpx    = zeros((3,nxzones,nyzones))
        self.qpy    = zeros((3,nxzones,nyzones))
        self.qmx    = zeros((3,nxzones,nyzones))
        self.qmy    = zeros((3,nxzones,nyzones))
        self.nx     = nxzones
        self.ny     = nyzones
        self.gx     = 3
        self.gy     = 3

    def setup_grid(self,xmin,xmax,ymin,ymax):
        dx = (xmax - xmin) / (self.nx - self.gx*2 - 1)
        dy = (ymax - ymin) / (self.ny - self.gy*2 - 1)
        xmin = xmin - self.gx*dx
        xmax = xmax + self.gx*dx
        ymin = ymin - self.gy*dy
        ymax = ymax + self.gy*dy
        for i in range(self.nx):
            # cell centers
            self.x[i] = xmin + (i)*dx
        for i in range(self.ny):
            self.y[i] = ymin + (i)*dy
        # cell LEFT interfaces
        for i in range(self.nx):
            self.xi[i] = self.x[i] - 0.5*dx
        for i in range(self.ny):
            self.yi[i] = self.y[i] - 0.5*dy

    
    def setup_ID(self):
        # Shocktube initial data
        rchange = (self.x[self.nx-self.gx-1] - self.x[self.gx]) / 2.0
        print rchange
        rho1 = 1.0
        rho2 = 0.1
        press1 = 1.0
        press2 = 0.125
        for i in range(self.nx):
            for j in range(self.ny):
                if self.x[i] < rchange:
                    self.rho[i,j] = rho1
                    self.press[i,j] = press1
                    self.eps[i,j] = press1 / (rho1) / (gamma - 1.0)
                    self.vel[i,j] = 0.0
                else:
                    self.rho[i,j] = rho2
                    self.press[i,j] = press2
                    self.eps[i,j] = press2 / (rho2) / (gamma - 1.0)
                    self.vel[i,j] = 0.0


##################################### some basic functions
def prim2con(rho,vel,eps):
    q = zeros((3,rho.shape[0],rho.shape[1]))
    q[0,:,:] = rho
    q[1,:,:] = rho * vel
    q[2,:,:] = rho * eps + \
        0.5 * rho * vel**2

    return q

def con2prim(q):
    rho = q[0,:,:]
    vel = q[1,:,:] / rho
    eps = q[2,:,:] / rho - 0.5*vel**2
    press = eos_press(rho,eps,gamma)

    return (rho,eps,press,vel)

############# boundary conditions
def apply_bcs(hyd):
    hyd.rho[0:hyd.gx-1] = hyd.rho[hyd.gx]
    hyd.vel[0:hyd.gx-1] = hyd.vel[hyd.gx]
    hyd.eps[0:hyd.gx-1] = hyd.eps[hyd.gx]
    hyd.press[0:hyd.gx-1] = hyd.press[hyd.gx]
    
    hyd.rho[0:hyd.gy-1] = hyd.rho[hyd.gy]
    hyd.vel[0:hyd.gy-1] = hyd.vel[hyd.gy]
    hyd.eps[0:hyd.gy-1] = hyd.eps[hyd.gy]
    hyd.press[0:hyd.gy-1] = hyd.press[hyd.gy]

    hyd.rho[hyd.nx-hyd.gx:hyd.nx-1] = hyd.rho[hyd.nx-hyd.gx-1]
    hyd.vel[hyd.nx-hyd.gx:hyd.nx-1] = hyd.vel[hyd.nx-hyd.gx-1]
    hyd.eps[hyd.nx-hyd.gx:hyd.nx-1] = hyd.eps[hyd.nx-hyd.gx-1]
    hyd.press[hyd.nx-hyd.gx:hyd.nx-1] = hyd.press[hyd.nx-hyd.gx-1]
    
    hyd.rho[hyd.ny-hyd.gy:hyd.ny-1] = hyd.rho[hyd.ny-hyd.gy-1]
    hyd.vel[hyd.ny-hyd.gy:hyd.ny-1] = hyd.vel[hyd.ny-hyd.gy-1]
    hyd.eps[hyd.ny-hyd.gy:hyd.ny-1] = hyd.eps[hyd.ny-hyd.gy-1]
    hyd.press[hyd.ny-hyd.gy:hyd.ny-1] = hyd.press[hyd.ny-hyd.gy-1]

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
        for i in range(hyd.gx-1,hyd.nx-hyd.gx+1):
            for j in range(hyd.gy-1,hyd.ny-hyd.gy+1):
                hyd.rhopx[i,j] = hyd.rho[i,j]
                hyd.rhomx[i,j] = hyd.rho[i,j]
                hyd.epspx[i,j] = hyd.eps[i,j]
                hyd.epsmx[i,j] = hyd.eps[i,j]
                hyd.velpx[i,j] = hyd.vel[i,j]
                hyd.velmx[i,j] = hyd.vel[i,j]
                hyd.rhopy[i,j] = hyd.rho[i,j]
                hyd.rhomy[i,j] = hyd.rho[i,j]
                hyd.epspy[i,j] = hyd.eps[i,j]
                hyd.epsmy[i,j] = hyd.eps[i,j]
                hyd.velpy[i,j] = hyd.vel[i,j]
                hyd.velmy[i,j] = hyd.vel[i,j]



    elif(type=='minmod'):
        for i in range(hyd.gy-1,hyd.ny-hyd.gy+1):
            (hyd.rhopx[:,i],hyd.rhomx[:,i]) = tvd_minmod_reconstruct(hyd.nx,hyd.gx,hyd.rho[:,i],hyd.x,hyd.xi)
            (hyd.epspx[:,i],hyd.epsmx[:,i]) = tvd_minmod_reconstruct(hyd.nx,hyd.gx,hyd.eps[:,i],hyd.x,hyd.xi)
            (hyd.velpx[:,i],hyd.velmx[:,i]) = tvd_minmod_reconstruct(hyd.nx,hyd.gx,hyd.vel[:,i],hyd.x,hyd.xi)
        for j in range(hyd.gx-1,hyd.nx-hyd.gx+1):
            (hyd.rhopy[j,:],hyd.rhomy[j,:]) = tvd_minmod_reconstruct(hyd.ny,hyd.gy,hyd.rho[j,:],hyd.y,hyd.yi)
            (hyd.epspy[j,:],hyd.epsmy[j,:]) = tvd_minmod_reconstruct(hyd.ny,hyd.gy,hyd.eps[j,:],hyd.y,hyd.yi)
            (hyd.velpy[j,:],hyd.velmy[j,:]) = tvd_minmod_reconstruct(hyd.ny,hyd.gy,hyd.vel[j,:],hyd.y,hyd.yi)

    elif(type=='mc'):
        for i in range(hyd.gy-1,hyd.ny-hyd.gy+1):
            (hyd.rhopx[:,i],hyd.rhomx[:,i]) = tvd_mc_reconstruct(hyd.nx,hyd.gx,hyd.rho[:,i],hyd.x,hyd.xi)
            (hyd.epspx[:,i],hyd.epsmx[:,i]) = tvd_mc_reconstruct(hyd.nx,hyd.gx,hyd.eps[:,i],hyd.x,hyd.xi)
            (hyd.velpx[:,i],hyd.velmx[:,i]) = tvd_mc_reconstruct(hyd.nx,hyd.gx,hyd.vel[:,i],hyd.x,hyd.xi)
        for j in range(hyd.gx-1,hyd.nx-hyd.gx+1):
            (hyd.rhopy[j,:],hyd.rhomy[j,:]) = tvd_mc_reconstruct(hyd.ny,hyd.gy,hyd.rho[j,:],hyd.y,hyd.yi)
            (hyd.epspy[j,:],hyd.epsmy[j,:]) = tvd_mc_reconstruct(hyd.ny,hyd.gy,hyd.eps[j,:],hyd.y,hyd.yi)
            (hyd.velpy[j,:],hyd.velmy[j,:]) = tvd_mc_reconstruct(hyd.ny,hyd.gy,hyd.vel[j,:],hyd.y,hyd.yi)
                
    else:
        print "reconstruction type not known; abort!"
        sys.exit()
    
    hyd.presspx = eos_press(hyd.rhopx,hyd.epspx,gamma)
    hyd.pressmx = eos_press(hyd.rhomx,hyd.epsmx,gamma)
    hyd.presspy = eos_press(hyd.rhopy,hyd.epspy,gamma)
    hyd.pressmy = eos_press(hyd.rhomy,hyd.epsmy,gamma)

    hyd.qpx = prim2con(hyd.rhopx,hyd.velpx,hyd.epspx)
    hyd.qmx = prim2con(hyd.rhomx,hyd.velmx,hyd.epsmx)
    hyd.qpy = prim2con(hyd.rhopy,hyd.velpy,hyd.epspy)
    hyd.qmy = prim2con(hyd.rhomy,hyd.velmy,hyd.epsmy)
    

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
    for i in range(hyd.gx,hyd.nx-hyd.gx):
        for j in range(hyd.gy,hyd.ny-hyd.gy):
            dtnew = min(dtnew, \
                        (hyd.x[i+1]-hyd.x[i]) / max(abs(hyd.vel[i,j]+cs[i,j]), \
                        abs(hyd.vel[i,j]-cs[i,j])), \
                        (hyd.y[j+1]-hyd.y[j]) / max(abs(hyd.vel[i,j]+cs[i,j]),\
                        abs(hyd.vel[i,j]-cs[i,j])))

    dtnew = min(cfl*dtnew,1.05*dtp)
    return dtnew

############# HLLE solver
def hlle(n,g,xi,rhop,rhom,velp,velm,qp,qm,pressp,pressm,epsp,epsm):
    fluxdiff = zeros((3,n))
    # compute eigenvalues
    evl  = zeros((3,n))
    evr  = zeros((3,n))
    smin = zeros(n)
    smax = zeros(n)
    csp  = sqrt(eos_cs2(rhop,epsp,gamma))
    csm  = sqrt(eos_cs2(rhom,epsm,gamma))
    for i in range(1,n-2):
        evl[0,i] = velp[i]
        evl[1,i] = velp[i] - csp[i]
        evl[2,i] = velp[i] + csp[i]
        evr[0,i] = velm[i+1]
        evr[1,i] = velm[i+1] - csm[i+1]
        evr[2,i] = velm[i+1] + csm[i+1]
        smin[i] = min(evl[0,i],evl[1,i],evl[2,i],evr[0,i],evr[1,i],evr[2,i],0.0)
        smax[i] = max(evl[0,i],evl[1,i],evl[2,i],evr[0,i],evr[1,i],evr[2,i],0.0)

    # set up flux left L and right R of the interface
    # at i+1/2
    fluxl = zeros((3,n))
    fluxr = zeros((3,n))

    
    for i in range(1,n-2):
        fluxl[0,i] = qp[0,i] * velp[i]
        fluxl[1,i] = qp[1,i] * velp[i] + pressp[i]
        fluxl[2,i] = (qp[2,i] + pressp[i]) * velp[i]
        fluxr[0,i] = qm[0,i+1] * velm[i+1]
        fluxr[1,i] = qm[1,i+1] * velm[i+1] + pressm[i+1]
        fluxr[2,i] = (qm[2,i+1] + pressm[i+1]) * velm[i+1]

    # solve the Riemann problem for the i+1/2 interface
    ds = smax - smin
    flux = zeros((3,n))
    for i in range(g-1,n-g+1):
        flux[:,i] = (smax[i] * fluxl[:,i] - smin[i] * fluxr[:,i] + smin[i] * smax[i] * (qm[:,i+1] - qp[:,i])) / (smax[i] - smin[i])

    # flux differences
    for i in range(g,n-g):
        dx = xi[i+1] - xi[i]
        fluxdiff[:,i] = 1.0/dx * (flux[:,i] - flux[:,i-1])
        
    return fluxdiff

############# HLLE solver
def hlle_2D(hyd):
    
    fluxdiffx = zeros((3,hyd.nx,hyd.ny))
    fluxdiffy = zeros((3,hyd.nx,hyd.ny))
    
    for i in range(hyd.gy-1,hyd.ny-hyd.gy+1):
        fluxdiffx[:,:,i] = hlle(hyd.nx,hyd.gx,hyd.xi,hyd.rhopx[:,i],hyd.rhomx[:,i],hyd.velpx[:,i],hyd.velmx[:,i],hyd.qpx[:,:,i],hyd.qmx[:,:,i],hyd.presspx[:,i],hyd.pressmx[:,i],hyd.epspx[:,i],hyd.epsmx[:,i])
        
    for j in range(hyd.gx-1,hyd.nx-hyd.gx+1):
        fluxdiffy[:,j,:] = hlle(hyd.ny,hyd.gy,hyd.yi,hyd.rhopy[j,:],hyd.rhomy[j,:],hyd.velpy[j,:],hyd.velmy[j,:],hyd.qpy[:,j,:],hyd.qmy[:,j,:],hyd.presspy[j,:],hyd.pressmy[j,:],hyd.epspy[j,:],hyd.epsmy[j,:])
        
    return fluxdiffx+fluxdiffy

############# RHS calculation
def calc_rhs(hyd):
    # reconstruction and prim2con
    hyd = reconstruct(hyd,reconstruction_type)
    # compute flux differences
    fluxdiff = hlle_2D(hyd)
    # return RHS = - fluxdiff
    return -fluxdiff

########################################################################
# Main program
########################################################################

# initialize
hyd = mydata(nxzones,nyzones)

# set up grid
hyd.setup_grid(0.0,1.0,0.0,0.5)

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
contourf(hyd.x,hyd.y,transpose(hyd.rho))
savefig("initial_data.pdf",format="pdf")

# main integration loop
while(t < tend):

    if(i % 10 == 0):
        # output
        print "%5d %15.6E %15.6E" % (i,t,dt)
        clf()
        contourf(hyd.x,hyd.y,transpose(hyd.rho))
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
contourf(hyd.x,hyd.y,transpose(hyd.rho))
xlabel('x')
ylabel('$\\rho$')
savefig("profile.pdf", format="pdf")

#f = open("mc.dat","w")
#for i in range(len(hyd.x)):
#    f.write("%15.6E %15.6E \n" % (hyd.x[i],hyd.rho[i]))
#f.close()
        


