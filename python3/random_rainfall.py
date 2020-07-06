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
# IMPORT PARAMETERS FROM CONFIGURATION FILE AND UNPACK
##################################################################

#spec = importlib.util.spec_from_file_location("config", sys.argv[1])
spec = importlib.util.spec_from_file_location("config","configs/config#2.py")
config = importlib.util.module_from_spec(spec)
spec.loader.exec_module(config)


## Rainfall
Rain0 = config.Rain0
rainfac = config.rainfac
rainpdf = config.rainpdf

##################################################################
# Random rainfall series in advance?
##################################################################
Nmeas = 50
trand1_series = np.random.uniform(0,1,Nmeas)
trand2_series = np.random.uniform(0,1,Nmeas)

# print('Uniform random numbers 1:', trand1_series)
# print('Uniform random numbers 2:', trand2_series)

nrain = np.empty(Nmeas)
Rr = np.empty(Nmeas)
Rm = np.empty(Nmeas)
nt = np.empty(Nmeas)

for i in range(0,Nmeas):
    if trand1_series[i] < 3/16:
        nrain[i] = 1 # Rain0
    elif trand1_series[i] < 10/16:
        nrain[i] = 2 # 2*Rain0
    elif trand1_series[i] < 15/16:
        nrain[i] = 4 # 4*Rain0
    else:
        nrain[i] = 9 # 9*Rain0

    ### "galton board 2": RAINFALL LOCATION
    if trand2_series[i] < 3/16:
        nloc = 1 # Reservoir
        Rm[i] = 0*Rain0 #moor
        Rr[i] = nrain[i]*Rain0 #res
        nt[i] = nrain[i]
    elif trand2_series[i] < 10/16:
        nloc = 2 # Moor and reservoir
        Rm[i] = nrain[i]*Rain0 #moor
        Rr[i] = nrain[i]*Rain0 #res
        nt[i] = 2*nrain[i]
    elif trand2_series[i] < 15/16:
        nloc = 3 # Moor
        Rm[i] = nrain[i]*Rain0 #moor
        Rr[i] = 0*Rain0 #res
        nt[i] = nrain[i]
    else:
        nloc = 4 # No rain in catchment
        Rm[i] = 0*Rain0 #moor
        Rr[i] = 0*Rain0 #res
        nt[i] = 0

print('Rain amounts', nrain)
print('Rm', Rm)
print('Rr', Rr)
print('nt', nt)

binedges = np.arange(0,20,1)

# the histogram of the data
# n, bins, patches = plt.hist(nt, bins=18, density=True, facecolor='b')
# print(bins)
hist, bin_edges = np.histogram(nt, bins=binedges, density=True)
print('hist', hist)
print('bins', bin_edges)
bin_edges = np.round(bin_edges,0)
plt.bar(bin_edges[:-1], hist, width = 1, color='#0504aa',alpha=0.7)
# plt.xlim(min(bin_edges), max(bin_edges))
plt.plot(rainfac,rainpdf,'ko')
plt.xlabel('Rainfall amount')
plt.ylabel('Probability')
plt.title('Histogram of rainfall amounts')
# plt.xlim(0, 18)
plt.show()


##################################################################
#
##################################################################
'''
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
    Rm = 0*Rain0*np.ones(1,Ny+1) #moor
    Rr = nrain*Rain0*np.ones(1,Ny+1) #res
    nt = nrain
elif trand2 < 10/16:
    nloc = 2 # Moor and reservoir
    Rm = nrain*Rain0*np.ones(1,Ny+1) #moor
    Rr = nrain*Rain0*np.ones(1,Ny+1) #res
    nt = 2*nrain
elif trand2 < 15/16:
    nloc = 3 # Moor
    Rm = nrain*Rain0*np.ones(1,Ny+1) #moor
    Rr = 0*Rain0*np.ones(1,Ny+1)#res
    nt = nrain
else:
    nloc = 4 # No rain in catchment
    Rm = 0*Rain0*np.ones(1,Ny+1) #moor
    Rr = 0*Rain0*np.ones(1,Ny+1) #res
    nt = 0
'''
