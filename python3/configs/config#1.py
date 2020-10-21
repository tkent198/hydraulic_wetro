#######################################################################
# Config file: parameters for Wetropolis Au dynamics
#######################################################################

'''
config#1: flood wave test.
Initialise with kinematic velocity and inflow and send through a Gaussian pulse.
'''

import numpy as np
from init_cond import init_cond_wetro0

## Output directory
outdir = '/configs/config#1'

## channel geometry: cross-section
hr = 0.015 #  depth river
wr = 0.05 # width
hf = 0.005 # flood plain depth
hc = 0.02 # depth city
hc = hr+hf # depth city
wf = 0.1  # width flood plain
wc = 0.1 # width city
tana = hf/wf # flood plain angle

## channel geometry: lengths and spatial coordinate
LR1 = 3.4 # length upto city region
LR2 = 3.8 # length from end of city region
LR3 = 4.2 # total length of channel
s_r = 0.932 #reservoir influx loc
s_m = 2.038 #moor influx loc
LR11 = 3.4 # transition zone from fp to c [LR11, LR1]
LR22 = 3.8 # transition zone from c to fp [LR2, LR22]
tr = 50 # severity of transition
Nk = int(25*LR3) # number of cells on computational mesh (25 times the domain length)

dbds = -0.01 # mean slope river bed

## other
g = 9.81     # acceleration of gravity
Cm = 0.02    # Manning coefficient
cfl = 0.5 # CFL number for stable time-stepping
Neq = 2 # U = (A, Au)
eta = 0.0004 #amplitude of Gaussian pulse for floodwave

# INITIAL and BOUNDARY CONDITIONS
# Periodic BC = 1
# Neumann BC = 2
# specified inflow BC = 3
ic = init_cond_wetro0
BC = 3

# time
tn = 0
wd = 10 #wetropolis day = 10s
tmax = 100
Nmeas = 100


## From WetropA1d.m:
# Spatial set-up of Wetropolis channel: length LR3 = circa 5m
# chan.L1 = 0.780;
# chan.r2 = 0.125;
# chan.L2 = pi*r2;
# chan.L3 = 0.637;
# chan.r4 = 0.175;
# chan.L4 = pi*r4;
# chan.L5 = 0.637;
# chan.r6 = 0.125;
# chan.L6 = pi*r6;
# chan.L7 = 0.427;
# chan.L10 = 0.094; # +del8;
# chan.L8 = 0.027;
# chan.del8 = L10-L8;
# chan.L7 = L7-del8;
# chan.L8 = L8+del8;
# chan.L9 = 0.463; # -del8;
# chan.r11 = 0.175;
# chan.L11 = 0.5*pi*r11;
# chan.L12 = 0.420;
# chan.LR1 = L1+L2+L3+L4+L5+L6+L7;
# chan.LR2 = LR1+L8+L9;
# chan.LR3 = LR2+L10+L11+L12;
