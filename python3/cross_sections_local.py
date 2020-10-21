##################################################################
# CHANNEL CROSS-SECTIONS ETC
##################################################################
import numpy as np


##################################################################

def xsec_Ahs(h,s,config):
    '''
    # function computes cross-section area A, wetted perimeter Wp, radius Rh
    # for channel geometry as a function of depth h and along-channel position s

    # Input:
    # (1) depth h
    # (2) coord s
    # (3) config

    # Output:
    # (1) area A
    # (2) wetted perimeter Wp
    # (3) hydraulic radius
    '''

    # unpack config pars
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

    if (s >= LR1) & (s <= LR2): # city region

        if (h < hc): # in rect channel
            area = h*wr
            Wp = wr + 2*h
        else: # > hc in flood
            area = h*(wr + 2*wc) - 2*wc*hc
            Wp = wr + 2*wc + 2*h


    elif (s > LR11) & (s < LR1): # transition from floodplain to city

        w = (s-LR11)/(LR1 - LR11) #linear
        w = 0.5*(1 + np.tanh(tr*(s - 0.5*(LR11+LR1)))) #smooth
        hrs = w*hc + (1-w)*hr
        hfs = hc - hrs
        tanas = hfs/wf

        if (h < hrs): # in rect channel
            area = h*wr
            Wp = wr + 2*h
        elif (h > hrs + hfs): #above slope
            area = h*(wr + wf) - wf*(hrs + 0.5*hfs)
            Wp = wr + 2*h - hfs + np.sqrt(hfs**2 + wf**2)
        else: # middle  sloped region
            area = h*wr + 0.5*(h - hrs)**2/tanas
            Wp = h + wr + hrs + (h - hrs)*np.sqrt(1 + tanas**-2)


    elif (s > LR2) & (s < LR22): # transition to floodplain from city

        w = (s-LR2)/(LR22 - LR2) #linear
        w = 0.5*(1 + np.tanh(tr*(s - 0.5*(LR11+LR1)))) #smooth
        hrs = w*hr + (1-w)*hc
        hfs = hc - hrs
        tanas = hfs/wf

        if (h < hrs): # in rect channel
            area = h*wr
            Wp = wr + 2*h
        elif (h > hrs + hfs): #above slope
            area = h*(wr + wf) - wf*(hrs + 0.5*hfs)
            Wp = wr + 2*h - hfs + np.sqrt(hfs**2 + wf**2)
        else: # middle  sloped region
            area = h*wr + 0.5*(h - hrs)**2/tanas
            Wp = h + wr + hrs + (h - hrs)*np.sqrt(1 + tanas**-2)


    else: # floodplain

        if (h < hr): # in rect channel
            area = h*wr
            Wp = wr + 2*h
        elif (h > hr + hf): #above slope
            area = h*(wr + wf) - wf*(hr + 0.5*hf)
            Wp = wr + 2*h - hf + np.sqrt(hf**2 + wf**2)
        else: # middle  sloped region
            area = h*wr + 0.5*(h - hr)**2/tana
            Wp = h + wr + hr + (h - hr)*np.sqrt(1 + tana**-2)

    Rh = area/Wp

    return area, Wp, Rh

##################################################################

def xsec_hAs(A,s,config):

    '''
    # function computes depth h, and derivative dh/dA
    # for channel geometry as a function of area A and along-channel position s

    # Input:
    # (1) area A
    # (2) coord s
    # (3) config file

    # Output:
    # (1) depth h
    # (2) derivative dhdA
    '''

    # unpack config pars
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


    # critical areas
    A1 = wr*hr   # threshold area river
    A2 = (hr+hf)*(wr+wf)-wf*(hr+0.5*hf) # 2nd threshold area river
    Ac = wr*hc # threshold area city

    if (s > LR1) & (s < LR2): # city region

        if (A < Ac): # in rect channel
            h = A/wr
            dhdA = 1/wr
        else: # > Ac in flood
            h = (A + 2*wc*hc)/(wr + 2*wc)
            dhdA = 1/(wr + 2*wc)


    elif (s > LR11) & (s < LR1): # transition from floodplain to city

        w = (s-LR11)/(LR1 - LR11)
        hrs = w*hc + (1-w)*hr
        hfs = hc - hrs
        tanas = hfs/wf

        A1 = wr*hrs   # threshold area river
        A2 = (hrs+hfs)*(wr+wf)-wf*(hrs+0.5*hfs) # 2nd threshold area river

        if (A < A1): # in rect channel
            h = A/wr
            dhdA = 1/wr
        elif (A > A2): #above slope
            h = (A + wf*(hrs + 0.5*hfs))/(wr + wf)
            dhdA = 1/(wr + wf)
        else: # middle  sloped region
            h = hrs - wr*tanas + np.sqrt(tanas**2*wr**2 + 2*(A - wr*hrs)*tanas)
            dhdA = tanas/np.sqrt(tanas**2*wr**2 + 2*(A - wr*hrs)*tanas)


    elif (s > LR2) & (s < LR22): # transition from city to floodplain

        w = (s-LR2)/(LR22 - LR2)
        hrs = w*hr + (1-w)*hc
        hfs = hc - hrs
        tanas = hfs/wf

        A1 = wr*hrs   # threshold area river
        A2 = (hrs+hfs)*(wr+wf) - wf*(hrs+0.5*hfs) # 2nd threshold area river

        if (A < A1): # in rect channel
            h = A/wr
            dhdA = 1/wr
        elif (A > A2): #above slope
            h = (A + wf*(hrs + 0.5*hfs))/(wr + wf)
            dhdA = 1/(wr + wf)
        else: # middle  sloped region
            h = hrs - wr*tanas + np.sqrt(tanas**2*wr**2 + 2*(A - wr*hrs)*tanas)
            dhdA = tanas/np.sqrt(tanas**2*wr**2 + 2*(A - wr*hrs)*tanas)


    else: # floodplain

        if (A < A1): # in rect channel
            h = A/wr
            dhdA = 1/wr
        elif (A > A2): #above slope
            h = (A + wf*(hr + 0.5*hf))/(wr + wf)
            dhdA = 1/(wr + wf)
        else: # middle  sloped region
            h = hr - wr*tana + np.sqrt(tana**2*wr**2 + 2*(A - wr*hr)*tana)
            # dhdA = sqrt(config.tana/(2*A)); WRONG!!!
            dhdA = tana/np.sqrt(tana**2*wr**2 + 2*(A - wr*hr)*tana)



    return h, dhdA

##################################################################
def plot_xsec_hAs(A,s,config):

    '''
    # function computes coords for plotting plots of water depth h
    # in cross section A at location s: h = h(A(s,t),s)

    # Input:
    # (1) area A
    # (2) coord s
    # (3) config parameters

    # Output: [X,Y,Xc,Yc,h]
    '''

    # unpack config pars
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

    # critical areas
    A1 = wr*hr   # threshold area river
    A2 = (hr+hf)*(wr+wf) - wf*(hr+0.5*hf) # 2nd threshold area river
    Ac = wr*hc # threshold area city

    if (s > LR1) & (s < LR2): # city region

        Xc = [-wc,-wc, 0, 0, wr, wr, wr+wc, wr+wc]
        Yc = [hc+hc, hc, hc, 0, 0, hc, hc, hc+hc]

        if (A < Ac): # in rect channel
            h = A/wr
            X = [0,0,wr,wr]
            Y = [h,0,0,h]
        else: # > Ac in flood
            h = (A + 2*wc*hc)/(wr + 2*wc)
            X = [-wc,-wc, 0, 0, wr, wr, wr+wc, wr+wc]
            Y = [h, hc, hc, 0, 0, hc, hc, h]

    else: # floodplain

        Xc = [0, 0, wr, wr, wr+wf, wr+wf]
        Yc = [hc+hc,0 ,0 ,hr, hr+hf, hc+hc]

        if (A < A1): # in rect channel
            h = A/wr
            X = [0,0,wr,wr]
            Y = [h,0,0,h]
        elif (A > A2): #above slope
            h = (A + wf*(hr + 0.5*hf))/(wr + wf)
            X = [0, 0, wr, wr, wr+wf, wr+wf]
            Y = [h,0 ,0 ,hr, hr+hf, h]
        else: # middle  sloped region
            h = hr - wr*tana + np.sqrt(tana**2*wr**2 + 2*(A - wr*hr)*tana)
            X = [0, 0, wr, wr, wr+(h-hr)/tana]
            Y = [h,0,0,hr,h]


    return X,Y,Xc,Yc,h
