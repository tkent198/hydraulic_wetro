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
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation as manimation
import pdb # Python DeBugging
import matplotlib.gridspec as gs
matplotlib.use("Agg")


##################################################################
# CUSTOM MODULES REQUIRED
##################################################################
from cross_sections_local import xsec_hAs, xsec_Ahs, plot_xsec_hAs

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

# # Check array shapes:
# print('Rr shape: ', Rr.shape)
# print('Rm shape: ', Rm.shape)
# print('U_array shape: ', U_array.shape)
# print('hcs shape: ', hcs.shape)
# print('timevec: ', timevec.shape)
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

fp = index_fp[60]
ct = index_city[5]


# with writer.saving(fig, "writer_test.mp4", 100):
#     for k in range(10):
#         # Create a new plot object
#         plt.plot(range(k), range(k), 'o')
#         writer.grab_frame()

# try:
#     writer = manimation.writers['avconv']
# except KeyError:
#     writer = manimation.writers['ffmpeg']

FFMpegWriter = manimation.writers['ffmpeg']
# metadata = dict(title='Movie Test', artist='Matplotlib', comment='Movie support!')
writer = FFMpegWriter(fps=5)#, metadata=metadata)

fig = plt.figure(figsize=(13, 6.5))


with writer.saving(fig, str(dirn+outdir+'mov.mp4'),100):

    while tn <= tmax:

        h =  h_array[:,:,T][0]
        A =  U_array[0,:,T]
        Au = U_array[1,:,T]

        # plt.ion() ##
        G = gs.GridSpec(3, 4)
        ax_rt = plt.subplot(G[0,0]) # Rain factor against t
        ax_rhist = plt.subplot(G[0,1]) # Rain hist
        ax_Aus = plt.subplot(G[0,2:]) # River: Au(s,t)
        ax_hm = plt.subplot(G[1,0]) # Moor: hm(y,t)
        ax_cr = plt.subplot(G[1,1]) # Canals and res.: h(t)
        ax_hs = plt.subplot(G[1, 2:]) # River: h(s,t)
        ax_Qh = plt.subplot(G[2,0]) # Rating curve: Q = Q(h)
        ax_hct = plt.subplot(G[2,1]) # h(s,t) in city with flood theshold
        ax_hxfp = plt.subplot(G[2,2]) # X-sec: fp
        ax_hxct = plt.subplot(G[2,3]) # X-sec: ct

        ## Rainfall: times series
        ax_rt.plot(timevec[:T+1],Rm[:T+1], marker = '$M$', linestyle = 'None')
        ax_rt.plot(timevec[:T+1],Rr[:T+1], marker = '$R$', linestyle = 'None')
        ax_rt.plot(timevec[:T+1],Rm[:T+1]+Rr[:T+1], marker = '$&$', linestyle = 'None')
        ax_rt.set_ylim(-0.5, 20)
        ax_rt.set_yticks(rainfac)
        ax_rt.set_yticklabels(rainfac)
        # axes[0,2].set_xlim(0, tmeasure)

        ## Rainfall: histogram
        hist, bin_edges = np.histogram(Rm[:T+1]+Rr[:T+1], bins = np.arange(0,20,1), density=True)
        # print('hist', hist)
        # print('bins', bin_edges)
        bin_edges = np.round(bin_edges,0)
        ax_rhist.bar(bin_edges[:-1], hist, width = 1, color='#0504aa',alpha=0.7)
        # plt.xlim(min(bin_edges), max(bin_edges))
        ax_rhist.plot(rainfac,rainpdf,'ko')
        ax_rhist.set_xlabel('Rainfall amount')
        ax_rhist.set_ylabel('Density')
        # axes[1,2].title('Histogram of rainfall amounts')
        ax_rhist.set_xlim(-1, 19)
        ax_rhist.set_xticks(rainfac)
        ax_rhist.set_xticklabels(rainfac)

        ## Moor
        ax_hm.plot(yy,hm_array[:,T])
        ax_hm.set_xlim([0,Ly])
        ax_hm.set_ylim([hmmin,hmmax])

        ## Q-h relationship in city (a la rating curve)
        # if (hct[1:]>hct[:-1]):
        #     axes[2,0].plot(hct,Qct,'2k')
        # else:
        #     axes[2,0].plot(hct,Qct,'1b')
        hct = h_array[0,ct+1,:T+1]
        Qct = U_array[1,ct+1,:T+1]
        ax_Qh.plot(hct[np.where(hct[1:]>hct[:-1])],Qct[np.where(hct[1:]>hct[:-1])],'2k')
        ax_Qh.plot(hct[np.where(hct[1:]<=hct[:-1])],Qct[np.where(hct[1:]<=hct[:-1])],'1b')
        ax_Qh.plot([hc,hc],[Qmin+0.0001,Qmax],'r:')
        ax_Qh.set_xlabel('h')
        ax_Qh.set_ylabel('Q')
        ax_Qh.set_xlim([hmin+0.01,hmax])
        ax_Qh.set_ylim([Qmin+0.0001,Qmax])
        ax_Qh.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        ax_Qh.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        ## Canals and res: time series
        ax_cr.plot(timevec[:T+1],hcs[0,:T+1], marker = '$1$', linestyle = 'None')
        ax_cr.plot(timevec[:T+1],hcs[1,:T+1], marker = '$2$', linestyle = 'None')
        ax_cr.plot(timevec[:T+1],hcs[2,:T+1], marker = '$3$', linestyle = 'None')
        ax_cr.plot(timevec[:T+1],h_res[:T+1]/10, marker = '$R$', linestyle = 'None')
        ax_cr.set_ylim([0.005,0.02])

        ## h(city) time series with flood threshold h_T
        ax_hct.plot([timevec[0], timevec[T]],[hc, hc],'r:')
        ax_hct.plot(timevec[:T+1],hct)
        ax_hct.set_ylim([0.01,0.03])

        ## h(s,t)
        ax_hs.plot([s[fp], s[fp]],[hmin,hmax],'r:')
        ax_hs.plot([s[ct], s[ct]],[hmin,hmax],'r:')
        ax_hs.fill([config.LR1, config.LR2,config.LR2,config.LR1], [0,0,config.hc,config.hc],'r',alpha=0.1,linestyle='None')
        ax_hs.plot([s_r, s_r],[hmin,hmax],'k:')
        ax_hs.plot([s_m, s_m],[hmin,hmax],'k:')
        ax_hs.plot([Lc1, Lc1],[hmin,hmax],'k:')
        ax_hs.set_ylim([hmin,hmax])
        ax_hs.set_xlim([0,L])
        ax_hs.set_ylabel('$h(s,t)$',fontsize=12)
        ax_hs.set_xlabel('$s$',fontsize=12)
        ax_hs.plot([s[:-1],s[1:]],[h_array[0,1:-1,T],h_array[0,1:-1,T]],'b', linewidth = 1.0)
        ax_hs.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        ## Au(s,t)
        ax_Aus.set_ylim([Qmin,Qmax])
        ax_Aus.set_xlim([0,L])
        ax_Aus.plot([s_r, s_r],[Qmin,Qmax],'k:')
        ax_Aus.plot([s_m, s_m],[Qmin,Qmax],'k:')
        ax_Aus.plot([Lc1, Lc1],[Qmin,Qmax],'k:')
        ax_Aus.set_ylabel('$Au(s,t)$',fontsize=12)
        ax_Aus.set_xlabel('$s$',fontsize=12)
        ax_Aus.plot([s[:-1],s[1:]],[U_array[1,1:-1,T],U_array[1,1:-1,T]],'b', linewidth = 1.0)
        ax_Aus.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

        ## h cross-section: floodplain
        X,Y,Xc,Yc,__ = plot_xsec_hAs(U_array[0,fp+1,T],s[fp],config)
        ax_hxfp.plot(Xc,Yc,'k', linewidth=2.0)
        ax_hxfp.fill(X,Y,'b',alpha=0.2)
        ax_hxfp.text(Xc[-1],0.5*config.hr,'$t=%.3g$' %timevec[T], fontsize=12, horizontalalignment='right')
        ax_hxfp.text(Xc[-1],0.25*config.hr,'$s=%.3g$' %s[fp],fontsize=12, horizontalalignment='right')

        ## h cross-section: city
        X,Y,Xc,Yc,__ = plot_xsec_hAs(U_array[0,ct+1,T],s[ct],config)
        ax_hxct.plot(Xc,Yc,'k', linewidth=2.0)
        ax_hxct.fill(X,Y,'b',alpha=0.2)
        ax_hxct.text(Xc[-1],0.5*config.hr,'$t=%.3g$' %timevec[T], fontsize=12, horizontalalignment='right')
        ax_hxct.text(Xc[-1],0.25*config.hr,'$s=%.3g$' %s[ct],fontsize=12, horizontalalignment='right')

        plt.tight_layout(pad=0.2, w_pad=0.01, h_pad=0.01)

        # plt.show(block=False)
        # plt.pause(0.01)
        # plt.close()

        writer.grab_frame()
        plt.cla()

        print(' ')
        print('Grab frame at T = : ', tn)
        print(' ')

        tn += dtmeasure
        T += 1
