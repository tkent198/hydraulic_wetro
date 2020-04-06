function [X,Y,Xc,Yc,h] = plot_xsec_hAs(A,s,geom,chan)

% function computes coords for plotting plots of water depth h 
% in cross section A at location s: h = h(A(s,t),s)

% Input: 
% (1) area A
% (2) coord s 
% (3) channel geometry pars  
% (4) channel parametrisation

% Output: [X,Y,Xc,Yc,h]


% critical areas
A1 = geom.wr*geom.hr;   % threshold area river
A2 = (geom.hr+geom.hf)*(geom.wr+geom.wf)-geom.wf*(geom.hr+0.5*geom.hf); % 2nd threshold area river
Ac = geom.wr*geom.hc; % threshold area city

if (s > chan.LR1) && (s < chan.LR2) % city region
    Xc = [-geom.wc,-geom.wc, 0, 0, geom.wr, geom.wr, geom.wr+geom.wc, geom.wr+geom.wc];
    Yc = [geom.hc+geom.hc, geom.hc, geom.hc, 0, 0, geom.hc, geom.hc, geom.hc+geom.hc];
    if (A < Ac) % in rect channel
        h = A/geom.wr;
        X = [0,0,geom.wr,geom.wr];
        Y = [h,0,0,h];
    else % > Ac in flood
        h = (A + 2*geom.wc*geom.hc)/(geom.wr + 2*geom.wc);
        X = [-geom.wc,-geom.wc, 0, 0, geom.wr, geom.wr, geom.wr+geom.wc, geom.wr+geom.wc];
        Y = [h, geom.hc, geom.hc, 0, 0, geom.hc, geom.hc, h];
    end
else % floodplain
    Xc = [0, 0, geom.wr, geom.wr, geom.wr+geom.wf, geom.wr+geom.wf];
    Yc = [geom.hc+geom.hc,0 ,0 ,geom.hr, geom.hr+geom.hf, geom.hc+geom.hc];
    if (A < A1) % in rect channel
        h = A/geom.wr;
        X = [0,0,geom.wr,geom.wr];
        Y = [h,0,0,h];
    elseif (A > A2) %above slope
        h = (A + geom.wf*(geom.hr + 0.5*geom.hf))/(geom.wr + geom.wf);
        X = [0, 0, geom.wr, geom.wr, geom.wr+geom.wf, geom.wr+geom.wf];
        Y = [h,0 ,0 ,geom.hr, geom.hr+geom.hf, h];
    else % middle  sloped region
        h = geom.hr - geom.wr*geom.tana + sqrt(geom.tana^2*geom.wr^2 + 2*(A - geom.wr*geom.hr)*geom.tana);
        X = [0, 0, geom.wr, geom.wr, geom.wr+(h-geom.hr)/geom.tana];
        Y = [h,0,0,geom.hr,h];
    end   
end

end