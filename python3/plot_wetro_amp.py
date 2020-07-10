#######################################################################
# Main plotting script for Wetropolis Au dynamics
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
from matplotlib import animation
import animatplot as amp # see ex_animatplot.py for example script.

##################################################################
# CUSTOM MODULES REQUIRED
##################################################################
from cross_sections import xsec_hAs, xsec_Ahs, plot_xsec_hAs

##################################################################
# IMPORT PARAMETERS FROM CONFIGURATION FILE
##################################################################
#spec = importlib.util.spec_from_file_location("config", sys.argv[1])
spec = importlib.util.spec_from_file_location("config","configs/config#2.py")
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
dirn = str(cwd+'/configs'+outdir)
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


fp = index_fp[50]
ct = index_city[5]

##################################################################
# Animation using animatplot
##################################################################

h =  h_array[:,:,:][0]
A =  U_array[0,:,:]
Au =  U_array[1,:,:]
fp = index_fp[50]
ct = index_city[5]


fig, axes = plt.subplots(2, 1, figsize=(10,6))

axes[0].plot([s[fp], s[fp]],[0,0.04],'k:')
axes[0].plot([s[ct], s[ct]],[0,0.04],'k:')
axes[0].fill([config.LR1, config.LR2,config.LR2,config.LR1],[0,0,config.hc,config.hc],'r',alpha=0.1,linestyle='None')
axes[0].set_ylim([0,0.03])
axes[0].set_xlim([0,L])
axes[0].set_ylabel('$h(s,t)$',fontsize=14)
# axes[0].set_xlabel('$s$',fontsize=14)

axes[1].set_ylim([0,0.0004])
axes[1].set_xlim([0,L])
axes[1].set_ylabel('$Au(s,t)$',fontsize=14)
axes[1].set_xlabel('$s$',fontsize=14)


hline_block = amp.blocks.Line(s,h[1:,:],ax=axes[0], t_axis=1, linewidth = 1.0, drawstyle='steps-post')

Auline_block = amp.blocks.Line(s,Au[1:,:],ax=axes[1], t_axis=1, linewidth = 1.0, drawstyle='steps-post')

timeline = amp.Timeline(timevec, fps=10)

# now to construct the animation
anim = amp.Animation([hline_block, Auline_block], timeline)
anim.controls()

# anim.save_gif('images/multiblock')
plt.show()
