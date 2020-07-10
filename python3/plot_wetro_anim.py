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


##################################################################
# Animation -- mostly working fine
##################################################################

# Problem: setting data to patch for cross-section plots
# Use a time-dep. line only in cross section.
# Strange issue (unresolved): setting data for piecewise-contsnat plots -- see plot_issue_setdata.py.

h =  h_array[:,:,:][0]
A =  U_array[0,:,:]
Au =  U_array[1,:,:]
fp = index_fp[50]
ct = index_city[5]

fig, axes = plt.subplots(2, 2, figsize=(10,6))

axes[0,0].plot([s[fp], s[fp]],[0,0.04],'k:')
axes[0,0].plot([s[ct], s[ct]],[0,0.04],'k:')
axes[0,0].fill([config.LR1, config.LR2,config.LR2,config.LR1],[0,0,config.hc,config.hc],'r',alpha=0.1,linestyle='None')
axes[0,0].set_ylim([0,0.04])
axes[0,0].set_xlim([0,L])
axes[0,0].set_ylabel('$h(s,t)$',fontsize=14)
axes[0,0].set_xlabel('$s$',fontsize=14)
h_line, = axes[0,0].plot([],[],'b', linewidth = 1.0, drawstyle='steps-post')
title = axes[0,0].set_title("")

axes[0,1].set_ylim([0,0.0006])
axes[0,1].set_xlim([0,L])
axes[0,1].set_ylabel('$Au(s,t)$',fontsize=14)
axes[0,1].set_xlabel('$s$',fontsize=14)
Au_line, = axes[0,1].plot([],[],'b', linewidth = 1.0, drawstyle='steps-post')

X,Y,Xc,Yc,__ = plot_xsec_hAs(A[fp+1,0],s[fp],config)
axes[1,0].plot(Xc,Yc,'k', linewidth=2.0)
axes[1,0].set_ylim([0,0.03])
# axes[1,0].set_xlim([0,L])
axes[1,0].set_ylabel('$h(s=%.2f,t)$' %s[fp],fontsize=14)
axes[1,0].set_xlabel('Cross-channel',fontsize=14)
hfp_fill, = axes[1,0].plot([],[],'b',linewidth=2.0)
# axes[1,0].text(Xc[-1],0.25*config.hr,'$s=%.3g$' %s[fp],fontsize=14, horizontalalignment='right')

X,Y,Xc,Yc,__ = plot_xsec_hAs(A[ct+1,0],s[ct],config)
axes[1,1].plot(Xc,Yc,'k', linewidth=2.0)
axes[1,1].set_ylim([0,0.03])
# axes[1,0].set_xlim([0,L])
axes[1,1].set_ylabel('$h(s=%.2f,t)$' %s[ct],fontsize=14)
axes[1,1].set_xlabel('Cross-channel',fontsize=14)
hct_fill, = axes[1,1].plot([],[],'b',linewidth=2.0)
# axes[1,1].text(0.15,0.25*config.hr,'$s=%.3g$' %s[ct],fontsize=14, horizontalalignment='right')


def init():
    h_line.set_data([], [])
    Au_line.set_data([], [])
    # hfp_line.set_data([], [])
    hfp_fill.set_data([], [])
    # hct_line.set_data([], [])
    hct_fill.set_data([], [])
    title.set_text(" ")
    # return (h_line, Au_line, hfp_line, hfp_fill, hct_line, hct_fill, title)
    return (h_line, Au_line, hfp_fill, hct_fill, title)


def animate(i):

    # for k in range(0,len(s)-1):
    #     h_line.set_data([s[k], s[k+1]],[h_array[0,k+1,i],h_array[0,k+1,i]])
    #     Au_line.set_data([s[k], s[k+1]],[Au[k+1,i],Au[k+1,i]])

    # h_line.set_data([s[:-1],s[1:]],[h[1:-1,i],h[1:-1,i]])
    h_line.set_data(s,h[1:,i])

    # Au_line.set_data([s[:-1],s[1:]],[Au[1:-1,i],Au[1:-1,i]])
    Au_line.set_data(s, Au[1:,i])

    X,Y,Xc,Yc,__ = plot_xsec_hAs(A[fp+1,i],s[fp],config)
    # hfp_line.set_data(Xc,Yc)
    hfp_fill.set_data([X[0], X[-1]],[Y[0], Y[-1]])


    X,Y,Xc,Yc,__ = plot_xsec_hAs(A[ct+1,i],s[ct],config)
    # hct_line.set_data(Xc,Yc)
    hct_fill.set_data([X[0], X[-1]],[Y[0], Y[-1]])

    title.set_text("t = %.2f" % timevec[i])

    # return (h_line, Au_line, hfp_line, hfp_fill, hct_line, hct_fill, title)
    return (h_line, Au_line, hfp_fill, hct_fill, title)

anim = animation.FuncAnimation(fig, animate, init_func=init,
    frames=len(timevec), interval=100, blit=True)

plt.show()

# print('s', len(s))
# print('h', len(h[:,0]))
# print('hline',h_line)
