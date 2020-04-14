#######################################################################
# Main run script for Wetropolis Au dynamics
#######################################################################

'''
Using matlab file AuNCP_wetro0.m as base
'''

##################################################################
# GENERIC MODULES REQUIRED
##################################################################
import numpy as np
import scipy as sp
import os
import errno
import sys
import importlib.util

##################################################################
# CUSTOM MODULES REQUIRED
##################################################################
from flux_function import NCPflux_Au
from cross_sections import xsec_hAs, xsec_Ahs

##################################################################
# IMPORT PARAMETERS FROM CONFIGURATION FILE
##################################################################
#spec = importlib.util.spec_from_file_location("config", sys.argv[1])
spec = importlib.util.spec_from_file_location("config","configs/config#1.py")
config = importlib.util.module_from_spec(spec)
spec.loader.exec_module(config)


## config pars
hr = config.hr
wr = config.wr
hf = config.hf
hc = hr+hf
wf = config.wf
wc = config.wc
tana = config.tana
LR1 = config.LR1
LR2 = config.LR2
LR3 = config.LR3
LR11 = config.LR11
LR22 = config.LR22
tr = config.tr
s_r = config.s_r
s_m = config.s_m
dbds = config.dbds
g = config.g
Cm = config.Cm
Neq = config.Neq
ic = config.ic
cfl = config.cfl
BC = config.BC

##################################################################
# Set up dirs
##################################################################
outdir = config.outdir
cwd = os.getcwd()
dirn = str(cwd+outdir)
try:
    os.makedirs(dirn)
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise

##################################################################
# Set up grid
##################################################################

L=LR3 #length of domain
Nk=25*L #number of gridcells (excluding ghost)
Nk = int(Nk)
Nf=Nk+1 #number of nodes
Kk=L/Nk #length of cell

s = np.linspace(0, L, Nk+1)
sBC = np.linspace(-Kk, L+Kk, Nk+3)  #node loc with ghosts

# locating floodplain/city
index_fp = np.where((s < LR1) | (s > LR2))
index_city = np.where((s >= LR1) & (s <= LR2))


##################################################################
# Initialise
##################################################################

U0, B, h0 = ic(s,Nk,config)

U0 = np.insert(U0,0,U0[:,0],axis=1)
U0 = np.insert(U0,-1,U0[:,-1],axis=1)
B = np.append(np.append(B[0], B),B[-1])
h0 = np.append(np.append(h0[0], h0),h0[-1])

Z0 = B + h0

##################################################################
# Define time parameters
##################################################################

tn = config.tn
wd = config.wd #wetropolis day
tmax = config.tmax
Nmeas = config.Nmeas
dtmeasure = tmax/Nmeas
tmeasure = dtmeasure
index = 1

##################################################################
# Define system arrays with ghost cells for BCs
##################################################################

# Arrays for integration
Flux = np.empty([Neq,Nk+1])
S = np.empty([Neq,Nk+2]) # for source terms (friction, bed slope, and rain)
UU = np.empty([Neq,Nk+2])
SL = np.empty(Nk+1)
SR = np.empty(Nk+1) # numerical speeds
VNC = np.empty([Neq,Nk+1]) # NCP integrals
area = np.empty(Nk+2)
Wp = np.empty(Nk+2)
Rh = np.empty(Nk+2)
dhdA = np.empty(Nk+2)

U = U0
h = h0

# Arrays for saving
U_array = np.empty((Neq,Nk+2,Nmeas+1))
U_array[:,:,0] = U0

h_array = np.empty((1,Nk+2,Nmeas+1))
h_array[:,:,0] = h0

Z_array = np.empty((1,Nk+2,Nmeas+1))
Z_array[:,:,0] = Z0

##################################################################
# Numerical integration from t0 to tmax
##################################################################

while tn < tmax:

    # numerical fluxes
    for j in range(0,Nk+1):
        Flux[:,j], SL[j], SR[j], VNC[:,j] = NCPflux_Au(U[:,j],U[:,j+1],s[j],config)

    #Determine hydraulic radius, h and dh/dA etc
    #ghosts
    h[0], dhdA[0] = xsec_hAs(U[0,0],0.5*(-Kk+0),config)
    h[-1], dhdA[-1] = xsec_hAs(U[0,-1],0.5*(L+L+Kk),config)
    area[0], Wp[0], Rh[0] = xsec_Ahs(h[0],0.5*(-Kk+0),config)
    area[-1], Wp[-1], Rh[-1] = xsec_Ahs(h[-1],0.5*(L+L+Kk),config)

    #interiors
    for j in range(1,Nk+1):
        h[j], dhdA[j] = xsec_hAs(U[0,j],0.5*(s[j-1]+s[j]),config)
        area[j], Wp[j], Rh[j] = xsec_Ahs(h[j],0.5*(s[j-1]+s[j]),config)

    # compute extraneous forcing terms S(U)
    S[0,:] = 0
    S[1,:] = -config.g*U[0,:]*config.dbds - config.g*config.Cm**2*U[1,:]*abs(U[1,:]/U[0,:])/Rh**(4/3)

    #determine timestep for stability using wave eigen-speeds
    lam1 = U[1,:]/U[0,:] + np.sqrt(config.g*U[0,:]*dhdA)
    lam2 = U[1,:]/U[0,:] - np.sqrt(config.g*U[0,:]*dhdA)
    maxlam = np.maximum(lam1,lam2)

    dt = config.cfl*min(Kk/maxlam)

    # update time given new time step
    tn = tn + dt

    if tn > tmeasure:
        dt = dt - (tn - tmeasure) + 1e-12
        tn = tmeasure + 1e-12

    # P fluxes as per the NCP theory
    Pp = 0.5*VNC + Flux
    Pm = -0.5*VNC + Flux

    # integrate forward to t+dt
    if (BC == 1): # periodic NOT UPDATED

        #UU = U - dt*(Pp(:,2:Nk+1) - Pm(:,1:Nk))./Kk + dt*S # NOTE: not updated
        print('Error: periodic BCs not programmed')

    elif (BC == 2): # neumann

        #interior
        UU[:,1:-1] = U[:,1:-1] - dt*(Pp[:,1:] - Pm[:,:-1])/Kk + dt*S[:,1:-1]
        #ghosts
        UU[:,0] = UU[:,1]
        UU[:,-1]= UU[:,-2]

    elif (BC == 3): # specified inflow

        # interior
        UU[:,1:-1] = U[:,1:-1] - dt*(Pp[:,1:] - Pm[:,:-1])/Kk + dt*S[:,1:-1]
        # ghosts
        UU[0,0] = U[0,1] # A -- is this OK??
        # UU(2,1) = U0(2,1)+0.00005*sin(tn/(4*pi)) % Au: sine wave
        UU[1,0] = U0[1,0] + 0.0004*np.exp(-((tn-0.25*tmax)**2)/50) # Au: exp pulse
        UU[:,-1] = UU[:,-2]

    # update arrays for A, Au and h
    U = UU

    #h ghosts
    h[0], __ = xsec_hAs(U[0,0],0.5*(-Kk+0),config)
    h[-1], __ = xsec_hAs(U[0,-1],0.5*(L+L+Kk),config)
    #h interior
    for j in range(1,Nk+1):
        h[j], __ = xsec_hAs(U[0,j],0.5*(s[j-1]+s[j]),config)

    if tn > tmeasure:

        U_array[:,:,index] = U
        h_array[:,:,index] = h
        Z_array[:,:,index] = h+B

        index = index + 1
        tmeasure = tmeasure + dtmeasure


##################################################################
# Reached tmax; save data and end.
##################################################################

print(' ')
print('***** DONE: end of simulation at time:', tn)
print(' ')

print(' Saving simulation data in:', dirn)

np.save(str(dirn+'/U_array'),U_array)
np.save(str(dirn+'/h_array'),h_array)
np.save(str(dirn+'/Z_array'),Z_array)
