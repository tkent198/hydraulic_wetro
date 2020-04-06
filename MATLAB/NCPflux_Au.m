function [Flux, SL, SR, VNC] = NCPflux_Au(UL,UR,s,geom,chan)

%%%----- Numerical NCP flux function -----%%%
%
% after Rhebergen et al. (2008)
% 
% Function gives the intercell NCP flux value for use in FV Godunov numerical
% scheme for the Wetropolis St Venant system:
%
% U_t + (F(U))_x + G(U) U_x = 0 with U = (A,Au), F = (Au, Au^2 + g h A)
% and G = [0 0; -gh 0];
%
% INPUT: 
% Ul = (Al,Aul) from cell k 
% Ur = (Ar,Aur) from cell k+1
% geom = geometry pars defined in AuNCP
% chan = channel pars defined in AuNCP
%
% OUTPUT:
% Flux vector following the NCP theory (Rhebergen et al. 2008), with
% associated numerical speeeds, and NCP integral term VNC

% work with (A,u) rather than (A,Au)
AL = UL(1); 
AR = UR(1);
if (AL < 1e-16)
    uL = 0;
else
    uL = UL(2)/AL;
end

if (AR < 1e-16)
    uR = 0;
else
    uR = UR(2)/AR;
end

% compute h = h(A,s) and its derivative wrt A
[hL, dhdAL] = xsec_hAs(AL,s-0.04,geom,chan);
[hR, dhdAR] = xsec_hAs(AR,s+0.04,geom,chan);

% compute left and right wave speeds from eigenvalues
SL = min(uL - sqrt(geom.g*AL*dhdAL), uR - sqrt(geom.g*AR*dhdAR));
SR = max(uL + sqrt(geom.g*AL*dhdAL), uR + sqrt(geom.g*AR*dhdAR));


% compute h = h(A,s) and for Gauss integral

% 2-point
zetapl = 0.5*(1+1/sqrt(3));
zetami = 0.5*(1-1/sqrt(3));

[h1, ~] = xsec_hAs(AL + zetami*(AR-AL),s,geom,chan);
[h2, ~] = xsec_hAs(AL + zetapl*(AR-AL),s,geom,chan);

VNC2 = -0.5*geom.g*(AR-AL)*(h1+h2);

% 3-point
% zetapl = 0.5*(1+sqrt(3/5));
% zetami = 0.5*(1-sqrt(3/5));
% 
% [h1, ~] = xsec_hAs(AL + zetami*(AR-AL),s,geom,chan);
% [h2, ~] = xsec_hAs(AL + zetapl*(AR-AL),s,geom,chan);
% [h0, ~] = xsec_hAs(AL + 0.5*(AR-AL),s,geom,chan);
% 
% VNC2 = -0.5*geom.g*(AR-AL)*((5/9)*h1 + (5/9)*h2 + (8/9)*h0);


% 7-point
% zeta=[-0.9491079123, -0.7415311855, -0.4058451513, 0.0, 0.4058451513, 0.7415311855, 0.9491079123];
% % Weighting coefficients
% w=[0.1294849661, 0.2797053914, 0.3818300505, 0.4179591836, 0.3818300505, 0.2797053914, 0.1294849661];
% 
% h = zeros(1,length(w));
% for i = 1:7
%     [h(i), ~] = xsec_hAs(AL + 0.5*(1+zeta(i))*(AR-AL),s,geom,chan);
% end
% 
% VNC2 = -0.5*geom.g*(AR-AL)*sum(w.*h);

VNC = [0; VNC2];  

% define flux depending on wave speed
if (SL > 0)
    FluxL = [AL*uL; AL*uL^2 + geom.g*hL*AL];
    Flux = FluxL - 0.5*VNC;
elseif (SR < 0)
    FluxR = [AR*uR; AR*uR^2 + geom.g*hR*AR];
    Flux = FluxR + 0.5*VNC;
elseif (SL < 0) && (SR > 0)
    FluxL = [AL*uL; AL*uL^2 + geom.g*hL*AL];
    FluxR = [AR*uR; AR*uR^2 + geom.g*hR*AR];
    FluxHLL = (FluxL*SR - FluxR*SL + SL*SR*(UR - UL))/(SR-SL);
    Flux = FluxHLL - (0.5*(SL+SR)/(SR-SL))*VNC;
else
    Flux = [0; 0];
end