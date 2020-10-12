
import numpy as np
from cross_sections import xsec_hAs, xsec_Ahs

def NCPflux_Au(UL,UR,s,config):

    '''
    ###----- Numerical NCP flux function -----###
    #
    # after Rhebergen et al. (2008)
    #
    # Function gives the intercell NCP flux value for use in FV Godunov numerical
    # scheme for the Wetropolis St Venant system:
    #
    # U_t + (F(U))_x + G(U) U_x = 0 with U = (A,Au), F = (Au, Au**2 + g h A)
    # and G = [0 0; -gh 0];
    #
    # Note: integrals computed with 2-pt Gauss;
    # 3- and 7-pt also coded below; comment out as desired.
    #
    # INPUT:
    # Ul = (Al,Aul) from cell k
    # Ur = (Ar,Aur) from cell k+1
    # config
    #
    # OUTPUT:
    # Flux vector following the NCP theory (Rhebergen et al. 2008), with
    # associated numerical speeeds, and NCP integral term VNC
    '''

#    Nk = len(s)-1
#    Kk = config.LR3/Nk
    Kk = 0.04

    # work with (A,u) rather than (A,Au)
    AL = UL[0]
    AR = UR[0]
    if (AL < 1e-16):
        uL = 0
    else:
        uL = UL[1]/AL

    if (AR < 1e-16):
        uR = 0
    else:
        uR = UR[1]/AR

    # compute h = h(A,s) and its derivative wrt A
    hL, dhdAL = xsec_hAs(AL,s-Kk,config)
    hR, dhdAR = xsec_hAs(AR,s+Kk,config)

    # compute left and right wave speeds from eigenvalues
    SL = min(uL - np.sqrt(config.g*AL*dhdAL), uR - np.sqrt(config.g*AR*dhdAR))
    SR = max(uL + np.sqrt(config.g*AL*dhdAL), uR + np.sqrt(config.g*AR*dhdAR))

    # compute h = h(A,s) Gauss integral
    # 2-point
    zetapl = 0.5*(1+1/np.sqrt(3))
    zetami = 0.5*(1-1/np.sqrt(3))

    h1, __ = xsec_hAs(AL + zetami*(AR-AL),s,config)
    h2, __ = xsec_hAs(AL + zetapl*(AR-AL),s,config)

    VNC2 = -0.5*config.g*(AR-AL)*(h1+h2)

    # 3-point
    # zetapl = 0.5*(1+np.sqrt(3/5))
    # zetami = 0.5*(1-np.sqrt(3/5))
    #
    # h1 = xsec_hAs(AL + zetami*(AR-AL),s,config)
    # h2 = xsec_hAs(AL + zetapl*(AR-AL),s,config)
    # h0 = xsec_hAs(AL + 0.5*(AR-AL),s,config)
    #
    # VNC2 = -0.5*config.g*(AR-AL)*((5/9)*h1 + (5/9)*h2 + (8/9)*h0)


    # 7-point
    # zeta=[-0.9491079123, -0.7415311855, -0.4058451513, 0.0, 0.4058451513, 0.7415311855, 0.9491079123];
    # # Weighting coefficients
    # w=[0.1294849661, 0.2797053914, 0.3818300505, 0.4179591836, 0.3818300505, 0.2797053914, 0.1294849661];
    #
    # h = np.empty(len(w))

    # for i in range(0,7):
    #     h[i]= xsec_hAs(AL + 0.5*(1+zeta[i])*(AR-AL),s,config)
    #
    #
    # VNC2 = -0.5*config.g*(AR-AL)*np.sum(w*h)

    VNC1 = 0
    VNC = np.array([VNC1, VNC2])

    # define flux depending on wave speed
    if (SL > 0):
        FluxL = np.array([AL*uL, AL*uL**2 + config.g*hL*AL])
        Flux = FluxL - 0.5*VNC
    elif (SR < 0):
        FluxR = np.array([AR*uR, AR*uR**2 + config.g*hR*AR])
        Flux = FluxR + 0.5*VNC;
    elif (SL < 0) and (SR > 0):
        FluxL = np.array([AL*uL, AL*uL**2 + config.g*hL*AL])
        FluxR = np.array([AR*uR, AR*uR**2 + config.g*hR*AR])
        FluxHLL = (FluxL*SR - FluxR*SL + SL*SR*(UR - UL))/(SR-SL)
        # # # Check HLL values for steady state case
        # print('HLL = ', FluxHLL)
        # print('F1L = ', AL*uL, '; F1R = ', AR*uR)
        # print('F2L = ', AL*uL**2 + config.g*hL*AL, '; F2R = ', AR*uR**2 + config.g*hR*AR)
        # print(' ')
        Flux = FluxHLL - (0.5*(SL+SR)/(SR-SL))*VNC
    else:
        Flux = np.zeros(len(UL))


    return Flux, SL, SR, VNC
