#######################################################################
# Plotting script for Wetropolis Au dynamics - test cases
#######################################################################

'''
Plots simulation data from main_wetro_tests.py

Specify in spec below:
> config#0 for steady state test
> config#1 for floodwave test

Output: non-interative plots that need closing manually to step forward in time. See plot_wetro_amp.py for better videos.
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
import matplotlib.pyplot as plt
from matplotlib import animation

##################################################################
# CUSTOM MODULES REQUIRED
##################################################################
from cross_sections import xsec_hAs, xsec_Ahs, plot_xsec_hAs

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
Nk = config.Nk
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
# Nk=25*L #number of gridcells (excluding ghost)
# Nk = int(Nk)
Nf=Nk+1 #number of nodes
Kk=L/Nk #length of cell

s = np.linspace(0, L, Nk+1)
sBC = np.linspace(-Kk, L+Kk, Nk+3)  #node loc with ghosts

# locating floodplain/city
index_fp = np.where((s < LR1) | (s > LR2))
index_fp = np.array(index_fp)[0]
index_city = np.where((s >= LR1) & (s <= LR2))
index_city = np.array(index_city)[0]


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


##################################################################
# Plotting at times T
##################################################################
T = tn

fp = index_fp[50]
ct = index_city[5]

# plt.ion() ## Note this correction

fig, axes = plt.subplots(2, 2, figsize=(12,8))

axes[0,0].plot([s[fp], s[fp]],[0,0.04],'k:')
axes[0,0].plot([s[ct], s[ct]],[0,0.04],'k:')
axes[0,0].fill([config.LR1, config.LR2,config.LR2,config.LR1],[0,0,config.hc,config.hc],'r',alpha=0.1,linestyle='None')
#text(0.85*xmax,0.9*hmax,['t=',num2str(tn)])
axes[0,0].set_ylim([0,0.04])
axes[0,0].set_xlim([0,L])
axes[0,0].set_ylabel('$h(s,t)$',fontsize=14)
axes[0,0].set_xlabel('$s$',fontsize=14)

axes[0,1].set_ylim([0,0.0006])
axes[0,1].set_xlim([0,L])
axes[0,1].set_ylabel('$Au(s,t)$',fontsize=14)
axes[0,1].set_xlabel('$s$',fontsize=14)

while T < tmax:

    h =  h_array[:,:,T][0]
    A =  U_array[0,:,T]
    Au =  U_array[1,:,T]

    # for k in range(0,len(s)-1):
    #     axes[0,0].plot([s[k], s[k+1]],[h[k+1],h[k+1]],'b', linewidth = 1.0)
    #     axes[0,1].plot([s[k], s[k+1]],[Au[k+1],Au[k+1]],'b', linewidth = 1.0)

    axes[0,0].plot([s[:-1],s[1:]],[h[1:-1],h[1:-1]],'b', linewidth = 1.0)
    axes[0,1].plot([s[:-1],s[1:]],[Au[1:-1],Au[1:-1]],'b', linewidth = 1.0)

    X,Y,Xc,Yc,__ = plot_xsec_hAs(A[fp+1],s[fp],config)
    axes[1,0].plot(Xc,Yc,'k', linewidth=2.0)
    axes[1,0].fill(X,Y,'b',alpha=0.1)
    axes[1,0].text(Xc[-1],0.5*config.hr,'$t=%.3g$' %timevec[T], fontsize=14, horizontalalignment='right')
    axes[1,0].text(Xc[-1],0.25*config.hr,'$s=%.3g$' %s[fp],fontsize=14, horizontalalignment='right')

    X,Y,Xc,Yc,__ = plot_xsec_hAs(A[ct+1],s[ct],config)
    axes[1,1].plot(Xc,Yc,'k', linewidth=2.0)
    axes[1,1].fill(X,Y,'b',alpha=0.1)
    axes[1,1].text(Xc[-1],0.5*config.hr,'$t=%.3g$' %timevec[T], fontsize=14, horizontalalignment='right')
    axes[1,1].text(Xc[-1],0.25*config.hr,'$s=%.3g$' %s[ct],fontsize=14, horizontalalignment='right')

    T += 1

    plt.show()
    plt.pause(0.1)
