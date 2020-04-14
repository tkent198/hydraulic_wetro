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
    if (s >= config.LR1) & (s <= config.LR2): # city region
        
        if (h < config.hc): # in rect channel
            area = h*config.wr
            Wp = config.wr + 2*h
        else: # > hc in flood
            area = h*(config.wr + 2*config.wc) - 2*config.wc*config.hc
            Wp = config.wr + 2*config.wc + 2*h
        
        
    elif (s > config.LR11) & (s < config.LR1): # transition from floodplain to city
        
        w = (s-config.LR11)/(config.LR1 - config.LR11) #linear
        w = 0.5*(1 + np.tanh(config.tr*(s - 0.5*(config.LR11+config.LR1)))) #smooth
        hrs = w*config.hc + (1-w)*config.hr
        hfs = config.hc - hrs
        tanas = hfs/config.wf
        
        if (h < hrs): # in rect channel
            area = h*config.wr
            Wp = config.wr + 2*h
        elif (h > hrs + hfs): #above slope
            area = h*(config.wr + config.wf) - config.wf*(hrs + 0.5*hfs)
            Wp = config.wr + 2*h - hfs + np.sqrt(hfs**2 + config.wf**2)
        else: # middle  sloped region
            area = h*config.wr + 0.5*(h - hrs)**2/tanas
            Wp = h + config.wr + hrs + (h - hrs)*np.sqrt(1 + tanas**-2)
          
        
    elif (s > config.LR2) & (s < config.LR22): # transition to floodplain from city
        
        w = (s-config.LR2)/(config.LR22 - config.LR2) #linear
        w = 0.5*(1 + np.tanh(config.tr*(s - 0.5*(config.LR11+config.LR1)))) #smooth
        hrs = w*config.hr + (1-w)*config.hc
        hfs = config.hc - hrs
        tanas = hfs/config.wf
        
        if (h < hrs): # in rect channel
            area = h*config.wr
            Wp = config.wr + 2*h
        elif (h > hrs + hfs): #above slope
            area = h*(config.wr + config.wf) - config.wf*(hrs + 0.5*hfs)
            Wp = config.wr + 2*h - hfs + np.sqrt(hfs**2 + config.wf**2)
        else: # middle  sloped region
            area = h*config.wr + 0.5*(h - hrs)**2/tanas
            Wp = h + config.wr + hrs + (h - hrs)*np.sqrt(1 + tanas**-2)
          
        
    else: # floodplain
        
        if (h < config.hr): # in rect channel
            area = h*config.wr
            Wp = config.wr + 2*h
        elif (h > config.hr + config.hf): #above slope
            area = h*(config.wr + config.wf) - config.wf*(config.hr + 0.5*config.hf)
            Wp = config.wr + 2*h - config.hf + np.sqrt(config.hf**2 + config.wf**2)
        else: # middle  sloped region
            area = h*config.wr + 0.5*(h - config.hr)**2/config.tana
            Wp = h + config.wr + config.hr + (h - config.hr)*np.sqrt(1 + config.tana**-2)
    
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
    
    # critical areas
    A1 = config.wr*config.hr   # threshold area river
    A2 = (config.hr+config.hf)*(config.wr+config.wf)-config.wf*(config.hr+0.5*config.hf) # 2nd threshold area river
    Ac = config.wr*config.hc # threshold area city
    
    if (s > config.LR1) & (s < config.LR2): # city region
        
        if (A < Ac): # in rect channel
            h = A/config.wr
            dhdA = 1/config.wr
        else: # > Ac in flood
            h = (A + 2*config.wc*config.hc)/(config.wr + 2*config.wc)
            dhdA = 1/(config.wr + 2*config.wc)
        
        
    elif (s > config.LR11) & (s < config.LR1): # transition from floodplain to city
        
        w = (s-config.LR11)/(config.LR1 - config.LR11)
        hrs = w*config.hc + (1-w)*config.hr
        hfs = config.hc - hrs
        tanas = hfs/config.wf
        
        A1 = config.wr*hrs   # threshold area river
        A2 = (hrs+hfs)*(config.wr+config.wf)-config.wf*(hrs+0.5*hfs) # 2nd threshold area river
        
        if (A < A1): # in rect channel
            h = A/config.wr
            dhdA = 1/config.wr
        elif (A > A2): #above slope
            h = (A + config.wf*(hrs + 0.5*hfs))/(config.wr + config.wf)
            dhdA = 1/(config.wr + config.wf)
        else: # middle  sloped region
            h = hrs - config.wr*tanas + np.sqrt(tanas**2*config.wr**2 + 2*(A - config.wr*hrs)*tanas)
            dhdA = tanas/np.sqrt(tanas**2*config.wr**2 + 2*(A - config.wr*hrs)*tanas)
           
        
    elif (s > config.LR2) & (s < config.LR22): # transition from city to floodplain
        
        w = (s-config.LR2)/(config.LR22 - config.LR2)
        hrs = w*config.hr + (1-w)*config.hc
        hfs = config.hc - hrs
        tanas = hfs/config.wf
        
        A1 = config.wr*hrs   # threshold area river
        A2 = (hrs+hfs)*(config.wr+config.wf)-config.wf*(hrs+0.5*hfs) # 2nd threshold area river
        
        if (A < A1): # in rect channel
            h = A/config.wr
            dhdA = 1/config.wr
        elif (A > A2): #above slope
            h = (A + config.wf*(hrs + 0.5*hfs))/(config.wr + config.wf)
            dhdA = 1/(config.wr + config.wf)
        else: # middle  sloped region
            h = hrs - config.wr*tanas + np.sqrt(tanas**2*config.wr**2 + 2*(A - config.wr*hrs)*tanas)
            dhdA = tanas/np.sqrt(tanas**2*config.wr**2 + 2*(A - config.wr*hrs)*tanas)
         
        
    else: # floodplain
        
        if (A < A1): # in rect channel
            h = A/config.wr
            dhdA = 1/config.wr
        elif (A > A2): #above slope
            h = (A + config.wf*(config.hr + 0.5*config.hf))/(config.wr + config.wf)
            dhdA = 1/(config.wr + config.wf)
        else: # middle  sloped region
            h = config.hr - config.wr*config.tana + np.sqrt(config.tana**2*config.wr**2 + 2*(A - config.wr*config.hr)*config.tana)
            # dhdA = sqrt(config.tana/(2*A)); WRONG!!!
            dhdA = config.tana/np.sqrt(config.tana**2*config.wr**2 + 2*(A - config.wr*config.hr)*config.tana)
           
        
    
    return h, dhdA