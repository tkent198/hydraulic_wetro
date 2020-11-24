#######################################################################
# Main run script for Wetropolis full system
#######################################################################

'''
Wetropolis v1. Numerical integration of fully-coupled components:
-- river channel (St Venant PDE)
-- canals (mass ODE)
-- reservoir (mass ODE)
-- moor (groundwater model PDE)
'''


SAVEDATA = 1 # 1 to save simulation data
LOADRAIN = 0 # 1 to load saved rain in config__.npz; 0 to generate random rainfall

##################################################################
# GENERIC MODULES REQUIRED
##################################################################
import numpy as np
import scipy as sp
import os
import errno
import sys
import pdb # Python DeBugging
from datetime import datetime
import importlib.util
import matplotlib.pyplot as plt
from matplotlib import animation

##################################################################
# CUSTOM MODULES REQUIRED
##################################################################
from flux_function import NCPflux_Au
from cross_sections_local import xsec_hAs, xsec_Ahs, plot_xsec_hAs

##################################################################
# IMPORT PARAMETERS FROM CONFIGURATION FILE AND UNPACK
# ---- see config file for parameter descriptions ----
##################################################################

#spec = importlib.util.spec_from_file_location("config", sys.argv[1])
spec = importlib.util.spec_from_file_location("config","configs/config#2.py")
config = importlib.util.module_from_spec(spec)
spec.loader.exec_module(config)


##################################################################
# Set up dirs
##################################################################
outdir = config.outdir
cwd = os.getcwd()
dirn = str(cwd+'/configs'+outdir)
try:
    os.makedirs(dirn)
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise

##################################################################
# Time
##################################################################

tn = config.tn
wd = config.wd #wetropolis day
tmax = config.tmax
Nmeas = config.Nmeas
dtmeasure = tmax/Nmeas
tmeasure = dtmeasure
index = 1

##################################################################
# Rainfall
##################################################################

Rain0 = config.Rain0
rainfac = config.rainfac
rainpdf = config.rainpdf

##################################################################
# River: St Venant system
##################################################################

##config parameters
g = config.g
cfl = config.cfl
Cf = config.Cf # constant in weir relations
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
Nk = config.Nk
s_r = config.s_r
s_m = config.s_m
dbds = config.dbds
Cm = config.Cm
Neq = config.Neq
ic = config.ic
BC = config.BC

## set up grid
L=LR3 #length of domain
# Nk=25*L #number of gridcells (excluding ghost)
# Nk = int(Nk)
Nf=Nk+1 #number of nodes
Kk=L/Nk #length of cell

s = np.linspace(0, L, Nk+1)
sBC = np.linspace(-Kk, L+Kk, Nk+3)  #node loc with ghosts

# locating floodplain/city regions
# locating floodplain/city
index_fp = np.where((s < LR1) | (s > LR2))
index_fp = np.array(index_fp)[0]
index_city = np.where((s >= LR1) & (s <= LR2))
index_city = np.array(index_city)[0]

fp = index_fp[60]
ct = index_city[5]

## Initial conditions
# using ic function specified in config
U0, B, h0 = ic(s,Nk,config)
# add ghost values
U0 = np.insert(U0,0,U0[:,0],axis=1)
U0 = np.insert(U0,-1,U0[:,-1],axis=1)
B = np.append(np.append(B[0], B),B[-1])
h0 = np.append(np.append(h0[0], h0),h0[-1])

Z0 = B + h0

## Define system arrays with ghost cells for BCs
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

# Arrays for saving
U_array = np.empty((Neq,Nk+2,Nmeas+1))
U_array[:,:,0] = U0

h_array = np.empty((1,Nk+2,Nmeas+1))
h_array[:,:,0] = h0

Z_array = np.empty((1,Nk+2,Nmeas+1))
Z_array[:,:,0] = Z0


##################################################################
# Moor: groundwater model (FD)
##################################################################
## config paramters
gam_m = config.gam_m # proportion of water entering canal from moor (zero in physical model)
Ly = config.Ly # length in m
wv = config.wv # width Hele-Shaw moor cell in m
sigm = config.sigm   # fraction of moor pores filled
sigm0 = config.sigm0  # fraction of water remaining in moor pores after water has seeped out
sigme = sigm  #
mpor = config.mpor   # porosity moor
nu = config.nu  # viscosity water m^2/s
kperm = config.kperm # permeability
alph = kperm/(mpor*nu*sigm)
### connecting (canal) channel between moor and river channel
hcm = config.hcm # initial depth
Pwm = config.Pwm # weir height
Lc = config.Lc # length

## set up grid
Ny = config.Ny #no. cells
dy = Ly/Ny #gridsize
yy = np.linspace(0, Ly, Ny+1) # grid
hm0 = 0*yy #initial cond

# array for saving
hm_array = np.empty((Ny+1,Nmeas+1))
hm_array[:,0] = hm0

##################################################################
# Canals: mass ODE
##################################################################
Lc3 = config.Lc3 # m distance to lock 3
Lc2 = config.Lc2  # m distance to lock 2
Lc1 = config.Lc1# distance along canal of first lock in m
Lsec3 = config.Lsec3 # length third canal section in m
Lsec2 = config.Lsec2 # length second canal section in m
Lsec1 = config.Lsec1  # length first canal section
wc1 = config.wc1        # width of canal in m
Pw3 = config.Pw3   # depth weir in canal section 3
Pw2 = config.Pw2     # depth weir in canal section 2
Pw1 = config.Pw1     # depth weir in canal section 1
canalmaxdepth = config.canalmaxdepth
Hcc3 = config.Hcc3 # dike height along canal segment 3
Hcc2 = config.Hcc2 # dike height along canal segment 2
Hcc1 = config.Hcc1 # dike height along canal segment 1
## initial depths
h1c0 = Pw1
h2c0 = Pw2
h3c0 = Pw3

##################################################################
# Reservoir: mass ODE
##################################################################
gam_r = config.gam_r # proportion of water entering canal from res
wres = config.wres # m
Lres = config.Lres # length reservoir in m
Pwr = config.Pwr # weir height
hres0 = Pwr # initial depth

##################################################################
# Initialise full system
##################################################################
U = U0
h = h0
hm = hm0
hres = hres0
h1c = h1c0
h2c = h2c0
h3c = h3c0

## Moor, res, canal 1 inflow location to river channel
nxsr = int(np.floor(s_r/Kk)) # river gridcell in which res water is added
nxsm = int(np.floor(s_m/Kk)) # river gridcell in which moor water is added
nxLc1 = int(np.floor(Lc1/Kk)) #river gridcell in which canal-1 water is added

# # Random rainfall series in advance? see random_rainfall.py for some tests
# trand1_series = np.random.uniform(0,1,Nmeas)
# trand2_series = np.random.uniform(0,1,Nmeas)

if (LOADRAIN == 1):
    print(' Using rainfall data from:', dirn)

    dat = np.load(str(dirn+outdir+'.npz'))
    Rr_load = dat['Rr']
    Rm_load = dat['Rm']

else:
    print(' Generating random rainfall data day-by-day.')

tunit = 0
ncc = 0

Rmfac = [] #Rain factor moor
Rrfac = [] #Rain factor res
Rtfac = [] #Rain factor total
Qct = U[1,ct+1] #discharge at city location
hct = h[ct+1] #water depth at city location
hc1 = h1c0 #water depth canal-1
hc2 = h2c0 #water depth canal-2
hc3 = h3c0 #water depth canal-3
h_res = hres0 #water depth Reservoir

##Plotting limits:
hmax = 0.03
hmin = 0
Qmax = 0.0004
Qmin = 0
hmmax = 0.12
hmmin = 0

##################################################################
# Numerical integration from t0 to tmax (forward Euler)
##################################################################

while tn < tmax:

    ##################################################################
    # RANDOM RAINFALL
    ##################################################################
    while tn >= tunit:

        print(' ')
        print('tunit: ', tunit)
        print(datetime.now())

        ncc += 1
        tunit += wd

        if (LOADRAIN == 0):

            print('Generating random rainfall... ')

            trand1 = np.random.uniform()
            trand2 = np.random.uniform()

            ### "galton board 1": RAINFALL AMOUNT
            if trand1 < 3/16:
                nrain = 1 # Rain0
            elif trand1 < 10/16:
                nrain = 2 # 2*Rain0
            elif trand1 < 15/16:
                nrain = 4 # 4*Rain0
            else:
                nrain = 9 # 9*Rain0

            ### "galton board 2": RAINFALL LOCATION
            if trand2 < 3/16:
                nloc = 1 # Reservoir
                Rm = 0*Rain0*np.ones(Ny+1) #moor
                Rr = nrain*Rain0*np.ones(Ny+1) #res
                nt = nrain
                print('Location: Reservoir only')
            elif trand2 < 10/16:
                nloc = 2 # Moor and reservoir
                Rm = nrain*Rain0*np.ones(Ny+1) #moor
                Rr = nrain*Rain0*np.ones(Ny+1) #res
                nt = 2*nrain
                print('Location: Moor and Reservoir')
            elif trand2 < 15/16:
                nloc = 3 # Moor
                Rm = nrain*Rain0*np.ones(Ny+1) #moor
                Rr = 0*Rain0*np.ones(Ny+1)#res
                nt = nrain
                print('Location: Moor only')
            else:
                nloc = 4 # No rain in catchment
                Rm = 0*Rain0*np.ones(Ny+1) #moor
                Rr = 0*Rain0*np.ones(Ny+1) #res
                nt = 0
                print('Location: none')
            print('Total factor: ', nt)

            Rmfac = np.append(Rmfac,Rm[0]/Rain0)
            Rrfac = np.append(Rrfac,Rr[0]/Rain0)
            Rtfac = np.append(Rtfac,nt)

        elif (LOADRAIN == 1):

            print('Loaded rainfall data... ')
            print(' ')
            Rm = Rain0*Rm_load[ncc-1]*np.ones(Ny+1)
            Rr = Rain0*Rr_load[ncc-1]*np.ones(Ny+1)
            print('Moor rainfall factor:        Rm = ', Rm_load[ncc-1])
            print('Reservoir rainfall factor:   Rr = ', Rr_load[ncc-1])

    ##################################################################
    # Time step: restrictions using wave speeds/velocities
    ##################################################################
    #River (St Venant): determine hydraulic radius, h and dh/dA etc
    #ghosts
    h[0], dhdA[0] = xsec_hAs(U[0,0],0.5*(-Kk+0),config)
    h[-1], dhdA[-1] = xsec_hAs(U[0,-1],0.5*(L+L+Kk),config)
    area[0], Wp[0], Rh[0] = xsec_Ahs(h[0],0.5*(-Kk+0),config)
    area[-1], Wp[-1], Rh[-1] = xsec_Ahs(h[-1],0.5*(L+L+Kk),config)

    #interiors
    for j in range(1,Nk+1):
        h[j], dhdA[j] = xsec_hAs(U[0,j],0.5*(s[j-1]+s[j]),config)
        area[j], Wp[j], Rh[j] = xsec_Ahs(h[j],0.5*(s[j-1]+s[j]),config)

    # River (kinematic)
    # Vr = (wr*hr/(2*hr+wr))**(2/3)*np.sqrt(-dbds)/Cm

    # Moor
    hmg = max(0.5*(hm[2:] + hm[1:-1] + 0.5*(hm[1:-1] + hm[:-2])))
    hmg = max(hmg,hm[-2]+hm[-1])

    # Time step: components
    # dtr = Kk/abs(Vr) # dt kinematic river
    dtm = 0.5*(dy**2)/(alph*g*max(hmg,0.001)) #dt moor
    lam1 = U[1,:]/U[0,:] + np.sqrt(g*U[0,:]*dhdA)
    lam2 = U[1,:]/U[0,:] - np.sqrt(g*U[0,:]*dhdA)
    maxlam = np.maximum(abs(lam1),abs(lam2))
    dtAu = min(Kk/maxlam) #dt st venant
    # dt = cfl*min(dtr,dtm,dtAu)# take min timestep with CFL
    dt = cfl*min(dtm,dtAu) # take min timestep with CFL

    ## debuggin time-step: resolved! Typo in config file for moor parameters
    # print('hmg: ', hmg)
    # print('dtm: ', dtm)
    # print('dtAu: ', dtAu)
    # print('dt: ', dt)
    # print(' ')
    #
    # pdb.set_trace()

    # update time given new time step
    tn = tn + dt

    if tn > tmeasure:
        dt = dt - (tn - tmeasure) + 1e-12
        tn = tmeasure + 1e-12

    ##################################################################
    # Reservoir
    ##################################################################
    hreso = hres
    Qresw = Cf*wres*np.sqrt(g)*max(hreso-Pwr,0)**(3/2)
    # Qresw = 0.0 # "infinite reservoir" -- sufficient for flood control?
    hres = hreso + dt*Rr[0] - (dt/(Lres*wres))*Qresw

    ##################################################################
    # Moor
    ##################################################################
    num = dt/(dy**2)
    hmo = hm

    ## first grid point
    # hm[0] = h3co #Canal-3 (?) Jan 2017: this is wrong needs correction h2co -> h3co -- why???
    # hm[0] = h[nxsm] #River at sm
    hm[0] = hcm # Moor canal

    ## interior points
    # hm(2:Ny) = hmo(2:Ny) + num*alph*g*(0.5*(hmo(3:Ny+1)+hmo(2:Ny)).*(hmo(3:Ny+1)-hmo(2:Ny))-0.5*(hmo(2:Ny)+hmo(1:Ny-1)).*(hmo(2:Ny)-hmo(1:Ny-1)))+dt*Rm(2:Ny)/(mpor*sigme);

    hm[1:-1] = hmo[1:-1] + num*alph*g*(0.5*(hmo[2:] + hmo[1:-1])*(hmo[2:] - hmo[1:-1]) - 0.5*(hmo[1:-1] + hmo[:-2])*(hmo[1:-1] - hmo[:-2])) + dt*Rm[1:-1]/(mpor*sigme)

    ## wall at last grid point so no flux
    hm[-1] = hmo[-1] + num*alph*g*(hmo[-2] + hmo[-1])*(hmo[-2] - hmo[-1]) + dt*Rm[-1]/(mpor*sigme)

    ## outflow
    Qmoor = 0.5*mpor*sigme*wv*alph*g*(hmo[1]**2 - hm[0]**2)/dy

    ## connecting canal
    hcmo = hcm
    Qcm = wr*np.sqrt(g)*Cf*max(hcmo-Pwm,0)**(3/2)
    hcm = hcmo + (dt/(Lc*wv))*( Qmoor - Qcm )

    ##################################################################
    # Canals
    ##################################################################
    h3co = h3c
    h2co = h2c
    h1co = h1c

    # Outflows Q3c, Q2c, Q1c
    Qc3 = wc1*np.sqrt(g)*Cf*max( h3co-Pw3, 0 )**(3/2)
    Qc2 = wc1*np.sqrt(g)*Cf*max( h2co-Pw2, 0 )**(3/2)
    Qc1 = wc1*np.sqrt(g)*Cf*max( h1co-Pw1, 0 )**(3/2)

    # first order approximations
    h3c = h3co + (dt/(Lsec3*wc1))*( gam_r*Qresw - Qc3 )
    h2c = h2co + (dt/(Lsec2*wc1))*( Qc3 + gam_m*Qmoor - Qc2 ) # note: gam_m = 0 in real set-up
    h1c = h1co + (dt/(Lsec1*wc1))*( Qc2-Qc1 )

    ##################################################################
    # River
    ##################################################################
    ## numerical fluxes
    for j in range(0,Nk+1):
        Flux[:,j], SL[j], SR[j], VNC[:,j] = NCPflux_Au(U[:,j],U[:,j+1],s[j],config)

    ## compute extraneous forcing terms S(U): bedslope and friction
    S[0,:] = 0
    S[1,:] = -g*U[0,:]*dbds - g*Cm**2*U[1,:]*abs(U[1,:]/U[0,:])/Rh**(4/3)

    ## P fluxes as per the NCP theory
    Pp = 0.5*VNC + Flux
    Pm = -0.5*VNC + Flux

    ## integrate forward to t+dt
    if (BC == 1): # periodic NOT UPDATED

        #UU = U - dt*(Pp(:,2:Nk+1) - Pm(:,1:Nk))./Kk + dt*S # NOTE: not updated
        print('Error: periodic BCs not programmed (nor suitable!)')

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
        UU[0,0] = U0[0,0] # A -- is this OK??
        # UU(2,1) = U0(2,1)+0.00005*sin(tn/(4*pi)) # Au: sine wave
        # UU[1,0] = U0[1,0] + 0.0004*np.exp(-((tn-0.25*tmax)**2)/50) # Au: exp pulse
        UU[1,0] = U0[1,0]
        UU[:,-1] = UU[:,-2]

    ## Area: update point source terms from moor, res, canal --- S_A
    UU[0,nxsm] = UU[0,nxsm] + dt*(1-gam_m)*Qcm/Kk # moor inflow
    UU[0,nxsr] = UU[0,nxsr] + dt*(1-gam_r)*Qresw/Kk # res inflow
    UU[0,nxLc1] = UU[0,nxLc1] + dt*Qc1/Kk # canal 1 inflow

    # Discharge: update point source term from moor, res, canal --- u*S_A
    UU[1,nxsm] = UU[1,nxsm] + (UU[1,nxsm]/UU[0,nxsm])*dt*(1-gam_m)*Qcm/Kk # moor inflow
    UU[1,nxsr] = UU[1,nxsr] + (UU[1,nxsr]/UU[0,nxsr])*dt*(1-gam_r)*Qresw/Kk # res inflow
    UU[1,nxLc1] = UU[1,nxLc1] + (UU[1,nxLc1]/UU[0,nxLc1])*dt*Qc1/Kk # canal 1 inflow

    # update arrays for A, Au and h
    U = UU

    #h ghosts
    h[0], __ = xsec_hAs(U[0,0],0.5*(-Kk+0),config)
    h[-1], __ = xsec_hAs(U[0,-1],0.5*(L+L+Kk),config)
    #h interior
    for j in range(1,Nk+1):
        h[j], __ = xsec_hAs(U[0,j],0.5*(s[j-1]+s[j]),config)

    if tn > tmeasure:

        # store data at fixed time intervals dtmeasure
        U_array[:,:,index] = U
        h_array[:,:,index] = h
        Z_array[:,:,index] = h+B
        hm_array[:,index] = hm
        Qct = np.append(Qct,U[1,ct+1])
        hct = np.append(hct,h[ct+1])
        hc1 = np.append(hc1,h1c)
        hc2 = np.append(hc2,h2c)
        hc3 = np.append(hc3,h3c)
        h_res = np.append(h_res, hres)

        index += 1
        tmeasure += dtmeasure

##################################################################
# Reached tmax; save data and end.
##################################################################

print(' ')
print('***** DONE: end of simulation at time:', tn)
print(' ')

if (LOADRAIN == 1):
    Rmfac = Rm_load
    Rrfac = Rr_load

if (SAVEDATA == 1):

    print('*** Saving simulation data in: ', dirn)

    np.savez(str(dirn+outdir), \
    U_array = U_array, \
    h_array = h_array, \
    Z_array = Z_array, \
    Rm = Rmfac,\
    Rr = Rrfac, \
    hcanals = np.vstack((hc1,hc2,hc3)), \
    h_res = h_res, \
    hm_array = hm_array)
else:
    print('*** Simulation data NOT saved ***')
