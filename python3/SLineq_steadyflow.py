#######################################################################
# Checking the wave speed inequalities -- non-trivial steady flow
#######################################################################


##################################################################
# GENERIC MODULES REQUIRED
##################################################################
import numpy as np
import scipy as sp
import os
import errno
import sys
import importlib.util
import matplotlib.pyplot as plt


##################################################################
# IMPORT PARAMETERS FROM CONFIGURATION FILE
##################################################################
#spec = importlib.util.spec_from_file_location("config", sys.argv[1])
spec = importlib.util.spec_from_file_location("config","configs/config#0.py")
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



##################################################################
# Define time parameters
##################################################################

tn = config.tn
wd = config.wd #wetropolis day
tmax = config.tmax
Nmeas = config.Nmeas
dtmeasure = tmax/Nmeas
tmeasure = dtmeasure
timevec = np.linspace(tn,tmax,Nmeas+1)

print(' Loading simulation data from:', dirn)

U_array = np.load(str(dirn+'/U_array.npy'))
h_array = np.load(str(dirn+'/h_array.npy'))
Z_array = np.load(str(dirn+'/Z_array.npy'))

h =  h_array[:,:,:][0]
A =  U_array[0,:,:]
Au =  U_array[1,:,:]

h =  h_array[0,0,0]
A =  U_array[0,0,0]
Au =  U_array[1,0,0]


Cm = 0.02
#hydraulic radius
RA = (A*wr)/(2*A + wr**2)
Rh = (h*wr)/(2*h + wr)

print(' A = ', A)
print(' h = ', h)
print(' RA = ', RA)
print(' Rh = ', Rh)

alpha1 = np.sqrt(-dbds*wr/g)/Cm
alpha2 = np.sqrt(-dbds/g)/Cm
print(' alpha1 = ', alpha1)
print(' alpha2 = ', alpha2)

# S^L,R as functions of a single fixed parameter
alpha1 = np.linspace(0, 2, 21)
alpha2 = np.linspace(0, 6, 21)
scSLA = alpha1*RA**(2/3) - np.sqrt(A)
scSLh = alpha2*Rh**(2/3) - np.sqrt(h)

fig, axes = plt.subplots(2, 1, figsize=(8,6))

axes[0].plot(alpha1,scSLA,'k')
axes[0].set_ylim(scSLA[0],scSLA[-1])
axes[0].set_xlim([alpha1[0],alpha1[-1]])
axes[0].set_ylabel('$\sqrt{wr/g} S^L$',fontsize=14)
axes[0].set_xlabel('$alpha1$',fontsize=14)
axes[0].axhline(y=0, color='k')

axes[1].plot(alpha2,scSLh,'k')
axes[1].set_ylim(scSLh[0],scSLh[-1])
axes[1].set_xlim([alpha2[0],alpha2[-1]])
axes[1].set_ylabel('$\sqrt{1/g} S^L$',fontsize=14)
axes[1].set_xlabel('$alpha2$',fontsize=14)
axes[1].axhline(y=0, color='k')

# S^L,R as functions of A and h?
h = np.linspace(0, 0.02, 21)
Rh = (h*wr)/(2*h + wr)
SLh = np.sqrt(-dbds)*Rh**(2/3)/Cm - np.sqrt(g*h)

A = h*wr
RA = (A*wr)/(2*A + wr**2)
SLA = np.sqrt(-dbds)*RA**(2/3)/Cm - np.sqrt(g*A/wr)
print('SLA = ', SLA)

fig, axes = plt.subplots(2, 1, figsize=(8,6))

axes[0].axhline(y=0, color='k')
axes[0].plot(A,SLA,'b')
axes[0].set_ylim(np.min(SLA),np.max(SLA))
axes[0].set_xlim([A[0],A[-1]])
axes[0].set_ylabel('$S^L$',fontsize=14)
axes[0].set_xlabel('$A$',fontsize=14)

axes[1].axhline(y=0, color='k')
axes[1].plot(h,SLh,'b')
axes[1].set_ylim(np.min(SLh),np.max(SLh))
axes[1].set_xlim([h[0],h[-1]])
axes[1].set_ylabel('$S^L$',fontsize=14)
axes[1].set_xlabel('$h$',fontsize=14)



# anim.save_gif('images/multiblock')
plt.show()
