function [area, Wp, Rh] = xsec_Ahs(h,s,geom,chan)

% function computes cross-section area A, wetted perimeter Wp, radius Rh
% for channel geometry as a function of depth h and along-channel position s

% Input: 
% (1) depth h
% (2) coord s 
% (3) channel geometry pars  
% (4) channel parametrisation

% Output: 
% (1) area A 
% (2) wetted perimeter Wp 
% (3) hydraulic radius

if (s > chan.LR1) && (s < chan.LR2) % city region
    
    if (h < geom.hc) % in rect channel
        area = h*geom.wr;
        Wp = geom.wr + 2*h;
    else % > hc in flood
        area = h*(geom.wr + 2*geom.wc) - 2*geom.wc*geom.hc;
        Wp = geom.wr + 2*geom.wc + 2*h;
    end
    
elseif (s > chan.LR11) && (s < chan.LR1) % transition from floodplain to city
    
    w = (s-chan.LR11)/(chan.LR1 - chan.LR11); %linear
    w = 0.5*(1 + tanh(chan.tr*(s - 0.5*(chan.LR11+chan.LR1)))); %smooth
    hrs = w*geom.hc + (1-w)*geom.hr;
    hfs = geom.hc - hrs;
    tanas = hfs/geom.wf;
    
    if (h < hrs) % in rect channel
        area = h*geom.wr;
        Wp = geom.wr + 2*h;
    elseif (h > hrs + hfs) %above slope
        area = h*(geom.wr + geom.wf) - geom.wf*(hrs + 0.5*hfs);
        Wp = geom.wr + 2*h - hfs + sqrt(hfs^2 + geom.wf^2);
    else % middle  sloped region
        area = h*geom.wr + 0.5*(h - hrs)^2/tanas;
        Wp = h + geom.wr + hrs + (h - hrs)*sqrt(1 + tanas^-2);
    end  
    
elseif (s > chan.LR2) && (s < chan.LR22) % transition to floodplain from city
    
    w = (s-chan.LR2)/(chan.LR22 - chan.LR2); %linear
    w = 0.5*(1 + tanh(chan.tr*(s - 0.5*(chan.LR11+chan.LR1)))); %smooth
    hrs = w*geom.hr + (1-w)*geom.hc;
    hfs = geom.hc - hrs;
    tanas = hfs/geom.wf;
    
    if (h < hrs) % in rect channel
        area = h*geom.wr;
        Wp = geom.wr + 2*h;
    elseif (h > hrs + hfs) %above slope
        area = h*(geom.wr + geom.wf) - geom.wf*(hrs + 0.5*hfs);
        Wp = geom.wr + 2*h - hfs + sqrt(hfs^2 + geom.wf^2);
    else % middle  sloped region
        area = h*geom.wr + 0.5*(h - hrs)^2/tanas;
        Wp = h + geom.wr + hrs + (h - hrs)*sqrt(1 + tanas^-2);
    end  
    
else % floodplain
    
    if (h < geom.hr) % in rect channel
        area = h*geom.wr;
        Wp = geom.wr + 2*h;
    elseif (h > geom.hr + geom.hf) %above slope
        area = h*(geom.wr + geom.wf) - geom.wf*(geom.hr + 0.5*geom.hf);
        Wp = geom.wr + 2*h - geom.hf + sqrt(geom.hf^2 + geom.wf^2);
    else % middle  sloped region
        area = h*geom.wr + 0.5*(h - geom.hr)^2/geom.tana;
        Wp = h + geom.wr + geom.hr + (h - geom.hr)*sqrt(1 + geom.tana^-2);
    end   
    
end

Rh = area/Wp;

end
        
    

