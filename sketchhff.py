####################################################################################
# 
####################################################################################
#
# GENERIC MODULES REQUIRED:
#
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from time import time
import os
import errno
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import time
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
from math import pi, e
import importlib.util
import sys
import math
import pandas as pd
from scipy.optimize import least_squares
#  from numpy import exp, sin, cosh, tanh, log, random
from lmfit import minimize, fit_report, Parameters
# 
plt.close("all")
#
# Definitions
#
def River_Param( xx, yy, nn1, nn2, sloc): # , data=None ): River_Param(xx,yy,n1,n2,ss)
    # River_Param is a function which parametrises the river using a local
    # coordinate, and then generates the global x, and y coordinates from this.
    # sloc (s-location) is a user input variable which corresponds to the
    # distance along the river.
    # The global coordinates x, and y are then computed from this local
    # parameter s, and are the output variables of this function.
    # The values given below correspond to the current physical model of
    # Wetropolis;
    #       Radii of the bends of the river, r
    #       Lengths of the straight sections of river, L
    #       The running total for distance along the river, s used to segregate
    #       the different sections of the river. 
    # The start of the river is taken to be the origin corresponding to
    # (x=0,y=0,s=0)
    # Accompanying this parametrisation code, is a figure titled
    # 'Wetropolis_river_param.png' which indicates how the river has been seperated into
    # sections.
    # The criteria used was differences in cross sectional area, and
    # geometrical properties.
    # Defining values corresponding to the current Wetropolis model
    # parameter characterizing the origin of the local coordinate frame
    x_s = 0.0
    #  radii of the curves
    r2 = 0.125
    r4 = 0.175
    r6 = 0.125 # same as r2
    r11 = 0.175 # same as r4
    # Array containing the lengths of the straight sections of the river
    # Rough estimates from Wout's routing data
    L = np.zeros(13)
    L[0] = x_s
    L[1] = 0.780
    L[3] = 0.637
    L[5] = 0.637
    L[7] = 0.427
    L[8] = 0.027
    L[9] = 0.463
    L[10] = 0.094
    L[12] = 0.420
    #  Array containing the s-values corresponding to the start/end of each
    #  section of river
    s = np.zeros(14)
    s[1] = x_s
    s[2] = s[1] + L[1]
    s[3] = s[2] + np.pi*r2
    s[4] = s[3] + L[3]
    s[5] = s[4] + np.pi*r4
    s[6] = s[5] + L[5]
    s[7] = s[6] + np.pi*r6
    s[8] = s[7] + L[7]
    s[9] = s[8] + L[8]
    s[10] = s[9] + L[9]
    s[11] = s[10] + L[10]
    s[12] = s[11] + 0.5*np.pi*r11
    s[13] = s[12] + L[12]
    #  Test the given value against conditions to find where on the river to look at
    #  condition loop
    #  loops over all of the values passed to the function in the variable sloc
    #  finds where the local coordinate is, and then generates global x and y
    #  coordinates
    print('sloc length=',0,len(sloc),s)
    #  plt.figure(16)
    for i in range(0,len(sloc)):
        if sloc[i] >= s[1] and sloc[i] < s[2]: #  first straight stretch of river - labelled as section (1)
            x = x_s + sloc[i]
            y = 0.0
            n1 = 0.0
            n2 = -1.0
        elif sloc[i] >= s[2] and sloc[i] <= s[3]: # first bend in the river - labelled as section (2)
            svec1 = sloc[i]
            x = x_s + L[1] +r2*np.sin((svec1 - s[2])/r2) # np.cos((svec1-s[2])/r2)
            y = r2*(1-np.cos((svec1-s[2])/r2)) # np.sin((svec1-s[2])/r2)
            n1 = np.sin((svec1 - s[2])/r2)
            n2 = -np.cos((svec1 - s[2])/r2)
            #  second straight stretch of river - labelled as section (3)
        elif sloc[i] > s[3] and sloc[i] <= s[4]:
            svec3 = sloc[i]
            x = x_s + L[1] -svec3 + s[3]
            y = 2*r2
            n1 = 0.0
            n2 = 1.0
        elif sloc[i] > s[4] and sloc[i] <= s[5]: #  second bend in the river - labelled as section (4)
            svec4 = sloc[i]
            x = x_s +L[1] - L[3] -r4*np.sin((svec4-s[4])/r4) # -np.cos
            y = 2*r2 + r4 - r4*np.cos((svec4-s[4])/r4) # np.sin
            n1 = np.sin((svec4-s[4])/r4)
            n2 = np.cos((svec4-s[4])/r4)
        elif sloc[i] > s[5] and sloc[i] <= s[6]: #  third straight stretch of river - labelled as section (5)
            svec5 = sloc[i]
            x = x_s +L[1]-L[3]+svec5-s[5]
            y = 2*(r2+r4)
            n1 = 0.0
            n2 = -1.0
        elif sloc[i] > s[6] and sloc[i] <= s[7]: #  third bend in the river - labelled as section (6)
            svec6 = sloc[i]
            x = x_s +L[1]-L[3]+L[5]+r6*np.sin((svec6-s[6])/r6) # np.cos((svec6-s[6])/r6)
            y = 2*(r2+r4)+r6-r6*np.cos((svec6-s[6])/r6) # np.sin((svec6-s[6])/r6)
            n1 = np.sin((svec6-s[6])/r6) 
            n2 = -np.cos((svec6-s[6])/r6) 
        elif sloc[i] > s[7] and sloc[i] <= s[8]: #  fourth straight stretch of river - labelled as section (7)
            svec7 = sloc[i]
            x = x_s +L[1]-L[3]+L[5]-(svec7-s[7])
            y = 2*(r2+r4+r6)
            n1 = 0.0
            n2 = 1.0
            #  fifth straight stretch of river - labelled as section (8)
            #  small area immediately before the city
        elif sloc[i] > s[8] and sloc[i] <= s[9]:
            svec8 = sloc[i]
            x = x_s +L[1]-L[3]+L[5]-L[7]-(svec8-s[8])
            y = 2*(r2+r4+r6)
            n1 = 0.0
            n2 = 1.0
            #  sixth straight stretch of river - labelled as section (9)
            #  city area
        elif sloc[i] > s[9] and sloc[i] <= s[10]:
            svec9 = sloc[i]
            x = x_s +L[1]-L[3]+L[5]-L[7]-L[8]-(svec9-s[9])
            y = 2*(r2+r4+r6);
            n1 = 0.0
            n2 = 1.0
        elif sloc[i] > s[10] and sloc[i] <= s[11]:
            #  seventh straight stretch of river - labelled as section (10)
            #  small area immediately after the city
            svec10 = sloc[i]
            x = x_s +L[1]-L[3]+L[5]-L[7]-L[8]-L[9]-(svec10-s[10])
            y = 2*(r2+r4+r6)
            n1 = 0.0
            n2 = 1.0
        elif sloc[i] > s[11] and sloc[i] <= s[12]:
            #  fourth bend of river - labelled as section (11)
            #  bend following the city
            svec11 = sloc[i]
            x = x_s + L[1]-L[3]+L[5]-L[7]-L[8]-L[9]-L[10]-r11*np.sin((svec11-s[11])/r11) # -np.cos((svec11-s[11])/r11)
            y = 2*(r2+r4+r6)+r11*np.cos((svec11-s[11])/r11)-r11 # -np.sin((svec11-s[11])/r11)
            n1 = -np.sin((svec11-s[11])/r11)
            n2 = np.cos((svec11-s[11])/r11)
            print(' i',i,x,y)
        elif sloc[i] > s[12] and sloc[i] <= s[13]:
            #  eighth straight stretch of river - labelled as section (12)
            svec12 = sloc[i]
            x = x_s + L[1]-L[3]+L[5]-L[7]-L[8]-L[9]-L[10]-r11*np.sin((s[12]-s[11])/r11)
            y = 2*(r2+r4+r6)-r11-(svec12-s[12])
            n1 = -1.0
            n2 = 0.0
            
        #  plt.plot(x,y,'.',lw=2)
        xx[i] = x
        yy[i] = y
        nn1[i] = n1
        nn2[i] = n2
        
    # print(' In param')
#
#
#

hf1 = 0.02
hf0 = 0.02
Lr1 = 3.0
Lr2 = 4.0
wf0 = 0.1
wcity = 0.1
hr0 = 0.015
hrc = 0.02
L = 5.2
as0 = 50.0
Ns = 520 # 250
ds = L/(1.0*Ns)
ss = np.linspace(0.0, L, Ns)
xx = 0.0*ss
yy = 0.0*ss
n1 = 0.0*ss
n2 = 0.0*ss
Nh = 200
wr = 0.05
hf = 0.005
tana = hf/wf0


hf = 0.5*hf1*(2.0-np.tanh(as0*(ss-Lr1))-np.tanh(as0*(Lr2-ss)))-0.5*hf0*(np.tanh(as0*(Lr1-ss))+np.tanh(as0*(ss-Lr2)))
wf = 0.5*wf0*(2.0-np.tanh(as0*(ss-Lr1))-np.tanh(as0*(Lr2-ss)))-wcity*(np.tanh(as0*(Lr1-ss))+np.tanh(as0*(ss-Lr2)))
#
# Plot
#
plt.figure(11)
plt.subplot(211)
plt.plot(ss,hf,'-k',lw=2)
plt.xlabel('$s$',fontsize=18)
plt.ylabel('$h_f(s)$',fontsize=18)
plt.subplot(212)
plt.plot(ss,wf,'-k',lw=2)
plt.xlabel('$s$',fontsize=18)
plt.ylabel('$w_f(s)$',fontsize=18)


LR1 = 3.8
LR2 = 4.2
LR11 = 3.6
LR22 = 4.4
hc = hrc
hr = hr0
dhr = 0.0
dhc = 0.02
tr = 50

hrs = 0.0*ss
dhrs = 0.0*ss
hfs = hrs
wfs = hrs
for nn in range(0, Ns):
    s = ss[nn]
    if (s >= LR1) & (s <=LR2 ): # city region
        hrs[nn] = hc
    elif (s > LR11) & (s < LR1): # transition from floodplain to city
        w = (s-LR11)/(LR1 - LR11) #linear
        w = 0.5*(1 + np.tanh(tr*(s - 0.5*(LR11+LR1)))) #smooth
        hrs[nn] = w*hc + (1-w)*hr
    elif (s > LR2) & (s < LR22): # transition to floodplain from city
        w = (s-LR2)/(LR22 - LR2) #linear
        w = 0.5*(1 + np.tanh(tr*(s - 0.5*(LR22+LR2)))) #smooth
        hrs[nn] = w*hr + (1-w)*hc
    else: # floodplain
        hrs[nn] = hr
        
w1 = (ss-LR11)/(LR1 - LR11) #linear
w1 = 0.5*(1 + np.tanh(tr*(ss - 0.5*(LR1+LR1)))) #smooth
w2 = (ss-LR2)/(LR22 - LR2) #linear
w2 = 0.5*(1 + np.tanh(tr*(ss - 0.5*(LR2+LR2)))) #smooth
w3 = 0.5*(1 + np.tanh(tr*(ss - LR11))) #smooth
w4 = 0.5*(1 + np.tanh(tr*(ss - LR1))) #smooth
w5 = 0.5*(1 + np.tanh(tr*(ss - LR2))) #smooth
w6 = 0.5*(1 + np.tanh(tr*(ss - LR22))) #smooth
dhrs = (w3-w4+w5-w6)*dhc #smooth
hrss = w1*hc + (1-w1)*hr - w2*(hc-hr) + (1-w2)*(hr-hr) #smooth
hfs = hc-hrss
wfs = w1*2*wcity + (1-w1)*wf0 - w2*(2*wcity-wf0) + (1-w2)*(wf0-wf0)
hrss = hrss+dhrs

#
# A(h,s) yields 3D plot; show at selected locations in s
# for nn in range(0, Nh):
#
plt.figure(15) 

ax = plt.axes(projection='3d')
ax.plot3D(ss,hrss,0,'--',lw=1)  # 
ax.plot3D(ss,hrss+hfs,0,'--',lw=1)  #
# s = 0, 1.25, 3.91, 4.5
hmax = 0.04
dh = hmax/(1.0*Nh)
hh = np.linspace(0.0, hmax, Nh)
for nn in range(0, Ns,1):
    s = ss[nn]
    area = 0*hh
    Wp = 0*hh
    h1rs = hrss[nn]
    h2fs = hfs[nn]
    for nh in range(0, Nh):
        # floodplain
        h = hh[nh]
        if (h < h1rs): # in rect channel
            area[nh] = h*wr
            Wp[nh] = wr + 2*h
        elif (h > h1rs + h2fs): #above slope
            area[nh] = h*(wr + wfs[nn]) - wfs[nn]*(h1rs + 0.5*h2fs)
            Wp[nh] = wr + 2*h - h2fs + np.sqrt(h2fs**2 + wfs[nn]**2)
        else: # middle sloped region
            area[nh] = h*wr + 0.5*(h - h1rs)**2/tana
            Wp[nh] = h + wr + h1rs + (h - h1rs)*np.sqrt(1 + tana**-2)
    ax.plot3D(s+0*hh,hh,area*1000,'-k',lw=1)  #
    
ax.set_box_aspect((20,4,4))
ax.set_xlim3d(0,L)
ax.set_ylim3d(0,hmax)
ax.set_zlim3d(0,0.006)
# ax.set_xticks([-2, 0, 2])
ax.set_yticks([0,0.02,0.04])
ax.set_zticks([0,2,4,6])
ax.set_xlabel('$s$')
ax.set_ylabel('$h$')
ax.set_zlabel('$A(h,s)x 10^{-3}$')


# 
#river bed plot x,y,z per s a cross profile on an orthogonal line
#
zfac = 100.0
slope = -0.01
zslope = slope*ss 
hb = np.max(hrss+hfs)
plt.figure(16)
ax1 = plt.axes(projection='3d')
print(' Before param',0)
River_Param(xx,yy,n1,n2,ss)
ax1.plot3D(xx,yy,0.0,'-k',lw=1)
ax1.plot3D(xx,yy,(-hb+zslope)*zfac,'-k',lw=1)
ax1.plot3D(xx+wr*n1,yy+wr*n2,(-hb+zslope)*zfac,'-k',lw=1)
ax1.plot3D(xx+wr*n1,yy+wr*n2,(-hb+hrss+zslope)*zfac,'-k',lw=1)
ax1.plot3D(xx+(wr+wfs)*n1,yy+(wr+wfs)*n2,(-hb+hrss+hfs+zslope)*zfac,'-k',lw=1)
ax1.plot3D(xx+(wr+wfs)*n1,yy+(wr+wfs)*n2,0.0,'-k',lw=1)
ax1.set_xlabel(' $x$ (m)',fontsize=14)
ax1.set_ylabel(' $y$ (m)',fontsize=14)
ax1.set_zlabel(' $z$ (cm)',fontsize=14)
ax1.set_zticks([-5,0])
# ax1.set_axis('square')
ax1.set_box_aspect((1,1*1.3/1.5,0.05))
ax1.set_xlim3d(-0.5,1)
ax1.set_ylim3d(-0.2,1)
ax1.set_zlim3d(-0.06*zfac,0.0*zfac)
# ax1.set_axis((-0.5,1,-0.5,1))
for nn in range(0, Ns,1):
    #  ax1.plot3D([xx[nn],xx[nn]+0.05*n1[nn]],[yy[nn],yy[nn]+0.05*n2[nn]],[0.0,0.0],'-',lw=2)
    ax1.plot3D([xx[nn],xx[nn],xx[nn]+wr*n1[nn],xx[nn]+wr*n1[nn],xx[nn]+(wr+wfs[nn])*n1[nn],xx[nn]+(wr+wfs[nn])*n1[nn]],[yy[nn],yy[nn],yy[nn]+wr*n2[nn],yy[nn]+wr*n2[nn],yy[nn]+(wr+wfs[nn])*n2[nn],yy[nn]+(wr+wfs[nn])*n2[nn]],[(0.0)*zfac,(-hb+zslope[nn])*zfac,(-hb+zslope[nn])*zfac,(-hb+hrss[nn]+zslope[nn])*zfac,(-hb+hrss[nn]+hfs[nn]+zslope[nn])*zfac,0.0],'-k',lw=1)
    # ax1.plot3D([xx[nn],xx[nn],xx[nn]+wr*n1[nn],xx[nn]+wr*n1[nn]], [yy[nn],yy[nn]-hb,yy[nn]+wr*n2[nn],yy[nn]+wr*n2[nn]], [0.0,-hb,-hb,-hb+hrss[nn]],'-',lw=2)
        
plt.figure(12)     
plt.subplot(221)
plt.plot(ss,hrs,'-k',lw=2)
plt.xlabel('$s$',fontsize=18)
plt.ylabel('$h_r(s)$',fontsize=18)
plt.subplot(222)
plt.plot(ss,hrss,'-k',lw=2)
plt.xlabel('$s$',fontsize=18)
plt.ylabel('Full $h_r(s)$',fontsize=18)
plt.subplot(223)
plt.plot(ss,hc-hrs,'-k',lw=2)
plt.xlabel('$s$',fontsize=18)
plt.ylabel('$h_f(s)$',fontsize=18)
plt.subplot(224)
plt.plot(ss,hfs,'-k',lw=2)
plt.xlabel('$s$',fontsize=18)
plt.ylabel('Full $h_f(s)$',fontsize=18)

plt.figure(13)     
plt.subplot(221)
plt.plot(ss,wfs,'-k',lw=2)
plt.xlabel('$s$',fontsize=18)
plt.ylabel('$w_f(s)$',fontsize=18)
plt.subplot(222)
plt.plot(ss,hrss,'-k',lw=2)
plt.xlabel('$s$',fontsize=18)
plt.ylabel('$h_r(s)$',fontsize=18)
plt.subplot(223)
plt.plot(ss,hfs,'-k',lw=2)
plt.xlabel('$s$',fontsize=18)
plt.ylabel('$h_f(s)$',fontsize=18)
plt.subplot(224)
plt.plot(ss,dhrs,'-k',lw=2)
plt.xlabel('$s$',fontsize=18)
plt.ylabel('$dhrs$',fontsize=18)

plt.show(block=True)
plt.pause(0.001)
plt.gcf().clear()
plt.show(block=False)

print("Finished program!")
