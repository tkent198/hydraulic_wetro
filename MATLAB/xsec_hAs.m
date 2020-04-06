function [h, dhdA] = xsec_hAs(A,s,geom,chan)

% function computes depth h, and derivative dh/dA
% for channel geometry as a function of area A and along-channel position s

% Input: 
% (1) area A
% (2) coord s 
% (3) channel geometry pars  
% (4) channel parametrisation

% Output: 
% (1) area h 
% (2) derivative dhdA

% critical areas
A1 = geom.wr*geom.hr;   % threshold area river
A2 = (geom.hr+geom.hf)*(geom.wr+geom.wf)-geom.wf*(geom.hr+0.5*geom.hf); % 2nd threshold area river
Ac = geom.wr*geom.hc; % threshold area city

if (s > chan.LR1) && (s < chan.LR2) % city region
    
    if (A < Ac) % in rect channel
        h = A/geom.wr;
        dhdA = 1/geom.wr;
    else % > Ac in flood
        h = (A + 2*geom.wc*geom.hc)/(geom.wr + 2*geom.wc);
        dhdA = 1/(geom.wr + 2*geom.wc);
    end
    
elseif (s > chan.LR11) && (s < chan.LR1)% transition from floodplain to city
    
    w = (s-chan.LR11)/(chan.LR1 - chan.LR11);
    hrs = w*geom.hc + (1-w)*geom.hr;
    hfs = geom.hc - hrs;
    tanas = hfs/geom.wf;
    
    A1 = geom.wr*hrs;   % threshold area river
    A2 = (hrs+hfs)*(geom.wr+geom.wf)-geom.wf*(hrs+0.5*hfs); % 2nd threshold area river
    
    if (A < A1) % in rect channel
        h = A/geom.wr;
        dhdA = 1/geom.wr;
    elseif (A > A2) %above slope
        h = (A + geom.wf*(hrs + 0.5*hfs))/(geom.wr + geom.wf);
        dhdA = 1/(geom.wr + geom.wf);
    else % middle  sloped region
        h = hrs - geom.wr*tanas + sqrt(tanas^2*geom.wr^2 + 2*(A - geom.wr*hrs)*tanas);
        dhdA = tanas/sqrt(tanas^2*geom.wr^2 + 2*(A - geom.wr*hrs)*tanas);
    end   
    
elseif (s > chan.LR2) && (s < chan.LR22)% transition from city to floodplain
    
    w = (s-chan.LR2)/(chan.LR22 - chan.LR2);
    hrs = w*geom.hr + (1-w)*geom.hc;
    hfs = geom.hc - hrs;
    tanas = hfs/geom.wf;
    
    A1 = geom.wr*hrs;   % threshold area river
    A2 = (hrs+hfs)*(geom.wr+geom.wf)-geom.wf*(hrs+0.5*hfs); % 2nd threshold area river
    
    if (A < A1) % in rect channel
        h = A/geom.wr;
        dhdA = 1/geom.wr;
    elseif (A > A2) %above slope
        h = (A + geom.wf*(hrs + 0.5*hfs))/(geom.wr + geom.wf);
        dhdA = 1/(geom.wr + geom.wf);
    else % middle  sloped region
        h = hrs - geom.wr*tanas + sqrt(tanas^2*geom.wr^2 + 2*(A - geom.wr*hrs)*tanas);
        dhdA = tanas/sqrt(tanas^2*geom.wr^2 + 2*(A - geom.wr*hrs)*tanas);
    end 
    
else % floodplain
    
    if (A < A1) % in rect channel
        h = A/geom.wr;
        dhdA = 1/geom.wr;
    elseif (A > A2) %above slope
        h = (A + geom.wf*(geom.hr + 0.5*geom.hf))/(geom.wr + geom.wf);
        dhdA = 1/(geom.wr + geom.wf);
    else % middle  sloped region
        h = geom.hr - geom.wr*geom.tana + sqrt(geom.tana^2*geom.wr^2 + 2*(A - geom.wr*geom.hr)*geom.tana);
%         dhdA = sqrt(geom.tana/(2*A)); WRONG!!!
        dhdA = geom.tana/sqrt(geom.tana^2*geom.wr^2 + 2*(A - geom.wr*geom.hr)*geom.tana);
    end   
    
end

end