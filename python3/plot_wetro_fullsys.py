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
import pdb # Python DeBugging


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

Rain0 = config.Rain0
rainfac = config.rainfac
rainpdf = config.rainpdf

Lc1 = config.Lc1

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

## Moor, res, canal 1 inflow location to river channel
nxsr = int(np.floor(s_r/Kk)) # river gridcell in which res water is added
nxsm = int(np.floor(s_m/Kk)) # river gridcell in which moor water is added
nxLc1 = int(np.floor(Lc1/Kk)) #river gridcell in which canal-1 water is added


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

dat = np.load(str(dirn+outdir+'.npz'))
print(' ')
print('Data files: ', dat.files)
print(' ')

## UNPACK
U_array = dat['U_array']
h_array = dat['h_array']
Rr = dat['Rr']
Rm = dat['Rm']
h_res = dat['h_res']
hm_array = dat['hm_array']
hcs = dat['hcanals']

## moor grid
Ly = config.Ly
yy = np.linspace(0, Ly, hm_array.shape[0]) # grid

# pdb.set_trace()

##################################################################
# Plotting at times T
##################################################################
T = tn

##Plotting limits:
hmax = 0.03
hmin = 0
Qmax = 0.0004
Qmin = 0
hmmax = 0.12
hmmin = 0

fp = index_fp[50]
ct = index_city[5]

# plt.ion() ## Note this correction


while T < tmax:

    h =  h_array[:,:,T][0]
    A =  U_array[0,:,T]
    Au = U_array[1,:,T]

    # plt.ion() ## Note this correction

    fig, axes = plt.subplots(3, 4, figsize=(13,7))#, constrained_layout=True)

    ## Rainfall: times series
    axes[0,0].plot(Rm[:T], marker = '$M$', linestyle = 'None')
    axes[0,0].plot(Rr[:T], marker = '$R$', linestyle = 'None')
    axes[0,0].plot(Rm[:T]+Rr[:T], marker = '$&$', linestyle = 'None')
    axes[0,0].set_ylim(-0.5, 20)
    axes[0,0].set_yticks(rainfac)
    axes[0,0].set_yticklabels(rainfac)
    # axes[0,2].set_xlim(0, tmeasure)

    ## Rainfall: histogram
    hist, bin_edges = np.histogram(Rm[:T]+Rr[:T], bins = np.arange(0,20,1), density=True)
    # print('hist', hist)
    # print('bins', bin_edges)
    bin_edges = np.round(bin_edges,0)
    axes[0,1].bar(bin_edges[:-1], hist, width = 1, color='#0504aa',alpha=0.7)
    # plt.xlim(min(bin_edges), max(bin_edges))
    axes[0,1].plot(rainfac,rainpdf,'ko')
    axes[0,1].set_xlabel('Rainfall amount')
    axes[0,1].set_ylabel('Density')
    # axes[1,2].title('Histogram of rainfall amounts')
    axes[0,1].set_xlim(-1, 19)
    axes[0,1].set_xticks(rainfac)
    axes[0,1].set_xticklabels(rainfac)


    ## Moor
    axes[1,0].plot(yy,hm_array[:,T])
    axes[1,0].set_xlim([0,Ly])
    axes[1,0].set_ylim([hmmin,hmmax])

    ## h-Q relationship in city (a la rating curve)
    # if (hct[1:]>hct[:-1]):
    #     axes[2,0].plot(hct,Qct,'2k')
    # else:
    #     axes[2,0].plot(hct,Qct,'1b')
    hct = h_array[0,ct+1,:T]
    Qct = U_array[1,ct+1,:T]
    axes[2,0].plot(hct[np.where(hct[1:]>hct[:-1])],Qct[np.where(hct[1:]>hct[:-1])],'2k')
    axes[2,0].plot(hct[np.where(hct[1:]<=hct[:-1])],Qct[np.where(hct[1:]<=hct[:-1])],'1b')
    axes[2,0].plot([hc,hc],[Qmin+0.0001,Qmax],'r:')
    axes[2,0].set_xlabel('h')
    axes[2,0].set_ylabel('Q')
    axes[2,0].set_xlim([hmin+0.01,hmax])
    axes[2,0].set_ylim([Qmin+0.0001,Qmax])
    axes[2,0].ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    axes[2,0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    ## Canals and res: time series
    axes[1,1].plot(hcs[0,:T], marker = '$1$', linestyle = 'None')
    axes[1,1].plot(hcs[1,:T], marker = '$2$', linestyle = 'None')
    axes[1,1].plot(hcs[2,:T], marker = '$3$', linestyle = 'None')
    axes[1,1].plot(h_res[:T]/10, marker = '$R$', linestyle = 'None')
    axes[1,1].set_ylim([0.005,0.015])

    ## h(city) time series with flood threshold h_T
    # axes[2,1].plot([0,ncc-1],[hc, hc],'r:')
    axes[2,1].plot(hct)
    axes[2,1].set_ylim([0,0.04])

    ## h(s,t)
    axes[0,2].plot([s[fp], s[fp]],[hmin,hmax],'r:')
    axes[0,2].plot([s[ct], s[ct]],[hmin,hmax],'r:')
    axes[0,2].fill([config.LR1, config.LR2,config.LR2,config.LR1], [0,0,config.hc,config.hc],'r',alpha=0.1,linestyle='None')
    axes[0,2].plot([s_r, s_r],[hmin,hmax],'k:')
    axes[0,2].plot([s_m, s_m],[hmin,hmax],'k:')
    axes[0,2].plot([Lc1, Lc1],[hmin,hmax],'k:')
    axes[0,2].set_ylim([hmin,hmax])
    axes[0,2].set_xlim([0,L])
    axes[0,2].set_ylabel('$h(s,t)$',fontsize=12)
    axes[0,2].set_xlabel('$s$',fontsize=12)
    axes[0,2].plot([s[:-1],s[1:]],[h_array[0,1:-1,T],h_array[0,1:-1,T]],'b', linewidth = 1.0)
    axes[0,2].ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    ## Au(s,t)
    axes[1,2].set_ylim([Qmin,Qmax])
    axes[1,2].set_xlim([0,L])
    axes[1,2].plot([s_r, s_r],[Qmin,Qmax],'k:')
    axes[1,2].plot([s_m, s_m],[Qmin,Qmax],'k:')
    axes[1,2].plot([Lc1, Lc1],[Qmin,Qmax],'k:')
    axes[1,2].set_ylabel('$Au(s,t)$',fontsize=12)
    axes[1,2].set_xlabel('$s$',fontsize=12)
    axes[1,2].plot([s[:-1],s[1:]],[U_array[1,1:-1,T],U_array[1,1:-1,T]],'b', linewidth = 1.0)
    axes[1,2].ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    ## h cross-section: floodplain
    X,Y,Xc,Yc,__ = plot_xsec_hAs(U_array[0,fp+1,T],s[fp],config)
    axes[2,2].plot(Xc,Yc,'k', linewidth=2.0)
    axes[2,2].fill(X,Y,'b',alpha=0.2)
    axes[2,2].text(Xc[-1],0.5*config.hr,'$t=%.3g$' %tmeasure, fontsize=12, horizontalalignment='right')
    axes[2,2].text(Xc[-1],0.25*config.hr,'$s=%.3g$' %s[fp],fontsize=12, horizontalalignment='right')

    ## h cross-section: city
    X,Y,Xc,Yc,__ = plot_xsec_hAs(U_array[0,ct+1,T],s[ct],config)
    axes[2,3].plot(Xc,Yc,'k', linewidth=2.0)
    axes[2,3].fill(X,Y,'b',alpha=0.2)
    axes[2,3].text(Xc[-1],0.5*config.hr,'$t=%.3g$' %tmeasure, fontsize=12, horizontalalignment='right')
    axes[2,3].text(Xc[-1],0.25*config.hr,'$s=%.3g$' %s[ct],fontsize=12, horizontalalignment='right')

    plt.tight_layout(pad=0.2, w_pad=0.01, h_pad=0.01)

    plt.show()
    plt.pause(0.1)

    T += 1
