##################################################################
#----------------- Initial conditions for wetropolis -----------------
##################################################################

'''
Functions generate different initial conditions described below for wetropolis

INPUT ARGS:
# x: mesh coords
# Neq: number of equations (variables) - 4 w/o topography, 5 w/ topography
# Nk: no. of cells in mesh
# H0: reference (scaled) height 
# L: length of domain
# A: amplitude
# V: velocity scale

OUTPUT:
# U0: array of initial data, size (Neq,Nk)

##################################################################
DESCRIPTIONS:

Rotation, no topography:

<init_cond_wetro0>
'''

###############################################################
import numpy as np
from cross_sections import xsec_Ahs
###############################################################

def init_cond_wetro0(s,Nk,config):

    '''
    ### Initial conditions for wetro0 system
    
    # INPUT:
    # s = cell edge coordinates
    # Nk = no. of cells
    # geom  = config
    
    # OUTPUT:
    # U0: [Neq, Nk] array consisting of FV [A,Au] data
    # B: [1, Nk] array consisting of FV topography data
    # h0: [1, Nk] array consisting of FV h data
    '''
    
    ###---------- ICs for the Neq equations ----------------###
    
    ### for depth h: flat
    ic_h = 0.0135*np.ones(len(s))
    # ic_h = 0.5*(2*config.hr+config.hf)*np.ones(len(s))
    # ic_h = 1.2*(config.hr+config.hf)*np.ones(len(s))
    
    ### for area A:
    ic_A = np.empty(len(s))
    W = np.empty(len(s))
    R = np.empty(len(s))
    
    for i in range(0,len(s)):
        ic_A[i], W[i], R[i] = xsec_Ahs(ic_h[i],s[i],config)
    
    
    ### kinetic solution for u(s,0)
    ic_u = R**(2/3)*np.sqrt(-config.dbds)/config.Cm
    
    ### for b: linear with gradient dbds
    L = config.LR3
    bc = -L*config.dbds
    B = np.linspace(bc,0,Nk+1)
    
    # Define array and fill with FV (piecewise constant) initial data 
    U0 = np.empty((config.Neq,Nk))
    B = 0.5*(B[0:Nk] + B[1:Nk+1]) # b
    h0 = 0.5*(ic_h[0:Nk] + ic_h[1:Nk+1]) # h
    U0[0,:] = 0.5*(ic_A[0:Nk]+ic_A[1:Nk+1]) # A
    U0[1,:] = 0.5*(ic_A[0:Nk]*ic_u[0:Nk] + ic_A[1:Nk+1]*ic_u[1:Nk+1]) # Au
    
    return U0, B, h0



###############################################################

#def init_cond_topog_cos(x,Nk,Neq,H0,L,A,V):
#    # superposition of cosines
#    ic1 = H0*np.ones(len(x))
#    ic2=1/ic1 # for hu = 1:
#    ic3 = np.zeros(len(x))
#
#    k = 2*np.pi
#    xp = 0.1
#    waven = [2,4,6]
#    A = [0.2, 0.1, 0.2]
#
#    B = A[0]*(1+np.cos(k*(waven[0]*(x-xp)-0.5)))+ A[1]*(1+np.cos(k*(waven[1]*(x-xp)-0.5)))+    A[2]*(1+np.cos(k*(waven[2]*(x-xp)-0.5)))
#    B = 0.5*B
#    
#    index = np.where(B<=np.min(B)+1e-10)
#    index = index[0]
#    B[:index[0]] = 0
#    B[index[-1]:] = 0
#
#    U0 = np.zeros((Neq,Nk))
#    B = 0.5*(B[0:Nk] + B[1:Nk+1]); # b
#    U0[0,:] = np.maximum(0, 0.5*(ic1[0:Nk] + ic1[1:Nk+1]) - B) # h
#    U0[1,:] = 0.5*(ic1[0:Nk]*ic2[0:Nk] + ic1[1:Nk+1]*ic2[1:Nk+1]) # hu
#    U0[2,:] = 0.5*(ic1[0:Nk]*ic3[0:Nk] + ic1[1:Nk+1]*ic3[1:Nk+1]) # hr
#    
#    return U0, B

###############################################################
