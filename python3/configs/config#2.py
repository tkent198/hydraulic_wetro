#######################################################################
# Config file: parameters for Wetropolis full system
#######################################################################
'''
Wetropolis v1: fixed parameters
-- general
-- river channel
-- canals
-- reservoir
-- moor (groundwater model)
'''

import numpy as np
from init_cond import init_cond_wetro0

##########################################
## Output directory
##########################################
outdir = '/config#2'

##########################################
## General
##########################################
g = 9.81     # acceleration of gravity
cfl = 0.5    # CFL number for stable time-stepping
Cf = (2/3)**(3/2) # constant in weir relations

##########################################
## Time
##########################################
tn = 0
wd = 10 #wetropolis day = 10s
Nmeas = 100
tmax = 1000


##########################################
## River
##########################################
Neq = 2 # U = (A, Au) for river dynamics
### cross-section channel geometry
hr = 0.015 #  depth river
wr = 0.05 # width
hf = 0.005 # flood plain depth
# hc = 0.02 # depth city
hc = hr+hf # depth city
wf = 0.1  # width flood plain
wc = 0.1 # width city
tana = hf/wf # flood plain angle

### spatial coordinate and lengths
LR1 = 3.8 # length upto city region
LR2 = 4.2 # length from end of city region
LR3 = 5.2 # total length of channel
LR11 = 3.6 # transition zone from fp to c [LR11, LR1]
LR22 = 4.4 # transition zone from c to fp [LR2, LR22]
tr = 50 # severity of transition

### coupling locations
s_r = 0.932 #reservoir influx loc
s_m = 2.038 #moor influx loc

### model parameters
dbds = -0.01 # mean slope river bed
Cm = 0.02    # Manning coefficient

### INITIAL and BOUNDARY CONDITIONS
ic = init_cond_wetro0
# Periodic BC = 1
# Neumann BC = 2
# specified inflow BC = 3
BC = 3

##########################################
## Canals
##########################################
Lc3 = 1.724 # m distance to lock 3
Lc2 = 3.608 # m distance to lock 2
Lc1 = 3.858 # distance along canal of first lock in m
Lsec3 = Lc3 # length third canal section in m
Lsec2 = Lc2-Lc3 # length second canal section in m
Lsec1 = Lc1-Lc2   # length first canal section
wc1 = 0.02         # width of canal in m
Pw3 = 0.0125      # depth weir in canal section 3
Pw2 = 0.0125      # depth weir in canal section 2
Pw1 = 0.01        # depth weir in canal section 1
canalmaxdepth = 0.02
Hcc3 = 0.06  # dike height along canal segment 2, canal max depth 0.0175m before overflow in river
Hcc2 = 0.04  # dike height along canal segment 2, canal max depth 0.0175m before overflow in river
Hcc1 = 0.021 # dike height along canal segment 1, canal max depth 0.0175m before overflow in river

##########################################
## Reservoir
##########################################
gam_r = 0.2 # proportion of water entering canal from reservoir (unknown)
### Geometry
wres = 0.123  # m
Lres  = 0.293    # length reservoir in m
Pwr  = 0.10 # weir height

##########################################
## Moor
##########################################
gam_m = 0.0 # proportion of water entering canal from moor (zero in physical model)
hr0 = 0.0135 # initial depth m?
### Geometry
Ly = 0.925 # length in m
wv = 0.095 # width Hele-Shaw moor cell in m

### Model parameters
sigm = 0.8    # fraction of moor pores filled
sigm0 = 0.1   # fraction of water remaining in moor pores after water has seeped out
sigme = sigm  #
mpor = 0.3    # porosity moor
nu = 10**(-6)  # viscosity water m^2/s
kperm = 10**(-8) # permeability
alph = kperm/(mpor*nu*sigm)

### connecting intermdeiate (canal) channel between moor and river channel
hcm = 0.0 # initial depth
Pwm = 0.02 # weir height
Lc = 0.1 # length


##########################################
## Rainfall
##########################################
Rain0 = 1.5*0.00013656 # This is r0 Eq. 17 in HESS article
rainfac = np.array([0,1,2,4,8,9,18])
rainpdf = np.array([16,24,77,89,35,8,7])/256
