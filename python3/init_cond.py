##################################################################
#----------------- Initial conditions for wetropolis -----------------
##################################################################


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
