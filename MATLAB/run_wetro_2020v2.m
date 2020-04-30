%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       WETROPOLIS RAINFALL AND FLOOD DEMONSTRATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - RIVER MODEL: ST. VENANT EQ.
% - GROUNDWATER MODEL: as before
% - RESERVOIR MODEL: as before
% - CANALS: as before

% v2 builds on run_wetro_2020.m. Updates here:

% > clean up code: lots of redundant parts and repetition now.
% > Advanced groundwater model (FEM) with connecting channel.
% > Live plotting improvements - save vids? Legend workaround?
% > More user friendly print text -- GUI with rainfall outcome?
% > EnKF for river model?

clear;
nparaflag = 2; % 1 is Onno's original design case; 2 is Luke/Onno/real case of parameter choice'
VIDEO = 0;
%% Canal
% Canal: s=0,..., Lc3 (lock 3), s=Lc3, ..., Lc2 (lock 2), s=Lc2,..., Lc1 (lock 1 into river)
% (a) Simple kinetic model with variables: h1=h1(t) in 1st, h2=h2(t) in 2nd section canal section
% (b) Shallow water model 1D or 2D?
% Solid wall at s=0.
%
%% Parameters
g = 9.81;          % m/s^2
Lc3 = 1.724; % m distance to lock 3
Lc2 = 3.608; % m distance to lock 2
Lc1 = 3.858; % distance along canal of first lock in m
Lsec3 = Lc3; % length third canal section in m
Lsec2 = Lc2-Lc3; % length second canal section in m
Lsec1 = Lc1-Lc2;   % length first canal section
wc = 0.02;         % width of canal in m
Cf = (2/3)^(3/2);
Pw3 = 0.0125;      % depth weir in canal section 3
Pw2 = 0.0125;      % depth weir in canal section 2
Pw1 = 0.01;        % depth weir in canal section 1
canalmaxdepth = 0.02;
Hcc3 = 0.06; % 0.04; % 0.06;    % dike height along canal segment 2, canal max depth 0.0175m before overflow in river
Hcc2 = 0.04; % 0.04; % 0.06;    % dike height along canal segment 2, canal max depth 0.0175m before overflow in river
Hcc1 = 0.021; % 0.022; % 0.032; % dike height along canal segment 1, canal max depth 0.0175m before overflow in river


%% River channel
%
% (a) simple river slope approximation: dh/dt	+ d Q(h)/dx = f
% Q(h) = h (wr*h/(2*h+wr))^2/3 sqrt(-db/dx)/Cm (simple rectangular river profile)
% wr*h/(2*h+wr) is area/wetted perimeter with wr width of river and h its depth
% (b) St Venant system (see TK Overleaf)
%
%% (a) kinematic model
Lx = 4.211; % m
dmean = 0.042; % 0.01*nfac; % was 0.04        % mean height difference across river bed in m
dbdxmean = -dmean/Lx; % mean sloep river bed
Cm = 0.02;            % Manning coefficient, see Table in Munson and ...; this is for glass
wr = 0.05;            % width of square river bed
Nx = 100;
dx = Lx/Nx;
xx = 0:dx:Lx;
Kk = xx(2:Nx+1)-xx(1:Nx);
xxm  = 0.5*(xx(2:Nx+1)+xx(1:Nx));
hr = 0*xx(1:Nx); % hr river depth defined as finite volumes, so one less than no grid points
Ar = wr*hr; % cross-sectional area for ractangular case
Vr = 0*hr;
bx = dmean*(Lx-xxm)/Lx; % sloping bed topography

%and canal
Hcx = Hcc3*0.5*(1+sign(Lc3-xxm))+Hcc2*0.5*(1+sign(xxm-Lc3)).*(1+sign(Lc2-xxm))*0.5+Hcc1*0.5*(1+sign(xxm-Lc2));
Hcxb = Hcx-canalmaxdepth; % Canal bottom topography (3 steps)
%
%
%% (b) St Venant river moidelling
%%%----- channel cross-section geometry  parameters ----%%%
geom.hr = 0.015; 
geom.wr = 0.05; % m width main river channel
geom.hf = 0.02; % m
geom.hf = 0.005;
geom.hc = 0.02;
geom.hc = geom.hr+geom.hf;
geom.wf = 0.1;  % m width flood plain
geom.wc = 0.1; 
geom.tana = geom.hf/geom.wf;
%%%----- other parameters ----%%%
geom.g = 9.81;     % acceleration of gravity
geom.Cm = 0.02;    % Manning coefficient
geom.dbds = -0.01; % mean slope river bed

%%%----- channel length parameters ----%%%
chan.LR3 = 4.2; % total length of channel
chan.LR3 = Lx; % total length of channel

chan.LR1 = 3.7; % length upto city region
chan.LR2 = 4.0; % length from end of city region

chan.LR11 = 3.4;% transition zone from fp to c [LR11, LR1]
chan.LR22 = 4.1;% transition zone from c to fp [LR2, LR22]

chan.s_r = 0.932; % reservoir influx
chan.s_m = 2.038; % moor influx

chan.tr = 50; % severity of transition 

%
%%%----- # of equations in system -----%%%
Neq = 2; % U = (A, Au)

%%%----- Set up grid with curvilinear coord s -----%%%
L=chan.LR3; %length of domain
Nk=100; %number of gridcells (excluding ghost)
Nf=Nk+1; %number of nodes
Kk=L/Nk; %length of cell
s=0:Kk:L;%node location
sBC = [-Kk s L+Kk]; %node loc with ghosts
ssm  = 0.5*(s(2:Nx+1)+s(1:Nx));% volume centre points (plotting)

% locating floodplain/city
index_fp = [find(s < chan.LR1) find(s > chan.LR2)];
index_city = intersect(find(s > chan.LR1), find(s < chan.LR2));
index_m = find(s > chan.s_m); 
index_r = find(s > chan.s_r); 

% BOUNDARY CONDITIONS
% Periodic BC = 1
% Neumann BC = 2
% specified inflow BC = 3
BC = 3; 

% ----- Apply initial conditions ----- %
[U0, B, h0] = initial_cond_wetro(s,Neq,Nk,geom,chan);

%with ghosts
U0 = [U0(:,1) U0 U0(:,end)];
B = [B(1) B B(end)];
h0 = [h0(1) h0 h0(end)];

Z0 = B + h0;

% Define empty system arrays for allocation  -- with ghost cells for BCs
Flux = zeros(Neq,Nk+1);
S = zeros(Neq,Nk+2); % for source terms (friction, bed slope, and rain)
UU = zeros(Neq,Nk+2);
SL = zeros(1,Nk+1); SR = zeros(1,Nk+1); % numerical speeds
VNC = zeros(Neq,Nk+1); % NCP integrals
area = zeros(1,Nk+2);
Wp = zeros(1,Nk+2);
Rh = zeros(1,Nk+2);
dhdA = zeros(1,Nk+2);

% take ICs
U = U0;
h = h0;

%% Groundwater model (moor)
% Moor water inflow in canal at s=Lv<Lc2 gam-fraction and river x=Lv (1-gam) -fraction with gam=gam(t) <=1
%	Orthogonal to canal
%       Option a) canal bottom coincide with Hele-Shaw moor bottom such that hm(1) = hc1
% hm = hm(y,t) for y=0,...,Ly
%% Parameters
switch nparaflag
    case 1 % Onnos original settings
        hr0 = 0.01; % m
        Ly = 0.4;   % m
        wv = 0.01;  % width Hele-Shaw moor cell in m
        Lv = 0.932; % m; % inflow of water off the Hele-Shaw moor in m
    case 2 % Luke/real settings
        hr0 = 0.0135; % m
        Ly = 0.925; % m
        wv = 0.095; % width Hele-Shaw moor cell in m
        Lv = 2.038; % m; % inflow of water off the Hele-Shaw moor in m
end

Ny = 20;
dy = Ly/Ny;
yy = 0:dy:Ly;  % Regular FD or FV grid
hm = 0*yy;     % Define moor depth variable on grid
sigm = 0.8;    % fraction of moor pores filled
sigm0 = 0.1;   % fraction of water remaining in moor pores after water has seeped out
sigme = sigm;  %
mpor = 0.3;    % porosity moor
nu = 10^(-6);  % viscosity water m^2/s
kperm = 10^-8; % permeability
alph = kperm/(mpor*nu*sigm);
gam = 0.2;


%% Reservoir
% Reservoir Hele-Shaw cell with adjustable weir
%
switch nparaflag
    case 1
        Lrer = 2.038;  % m
        wres = 0.01;   % m	%
        Lyy  = 0.4;    % length reservoir in
    case 3
        Lrer = 2.038;  % m
        wres = 0.01;   % m	%
        Lyy  = 0.4;    % length reservoir in
    case 2
        Lrer = 0.932;  % m
        wres = 0.123;  % m	%
        Lyy  = 0.293;    % length reservoir in m
end
Pwr  = 0.10; % weir height
hres = 0; % initially empty res

%% Moor, res, canal 1 coupling location
% (v)  canal weirs at s=Lc1 and s=Lc2 critical flow
%
%
nrer = ceil(Lrer/dx); %  grid box in which res water is added in the river
nxLv = ceil(Lv/dx);    % grid box in which moor water is added in the river
nxLc1 = ceil(Lc1/dx);  % grid box for regular grid in which canal water segment one is added
%

%% Rainfall and time
switch nparaflag
    case 1 % Onnos original
        Rain0 = 0.00175; % m/s
        Rain0 = 0.003; % 0.002 works m/s
        Rain = Rain0;
    case 2
        Rain0 = 0.00013656; % 0.002 works m/s % too low
        Rain0 = 1.5*0.00013656; % This is r0 Eq. 17 in HESS article
        Rain = Rain0;
end
%
TimeUnit = 10; % wd % 9 hrs from Kildwick till Armley times 3 is 27hrs; 12-20s across
Train = Lc2*wc*0.01/(Ly*wv*Rain); % TK: what is Train? unused (in case 2 below)
%
%
Vrate = Ly*wv*Rain0*[1,2,4,8,18]*TimeUnit*1000; % flow rates in liters/TimeUnit (Eq. 18 in HESS)


Tend = 100*TimeUnit;
dtmeas = 1.*TimeUnit;
tmeas = dtmeas;
tijd = 0;
% nmeas = dtmeas/dt;
%
%
%% Initialisation
% - canals
% - moor
% - kinematic river

figure(11); clf;
ini = 1;
switch ini
    case 1 % Mass check starting from absolute drought conditions and then it starts raining
        % canal sections
        h1c = 0;
        h2c = 0;
        h3c = 0;
        % Hele-Shaw moor
        hm = 0*yy;
        % River
        Ar = hr0*wr*ones(1,Nx);
end
figure(14); clf;
hro = h(1);
hrLo = h(nxLc1);
tijdo = tijd;
%
%
%
%
% nc1 = 0;
% nc2 = 0;
% nc4 = 0;
% nc8 = 0;
% nc16 = 0;
nt = 0;
nswitch = 0;
nswitch2 = 1;
tunit = 1;
nrain = 4;
Rm(1:Ny+1) = nrain*Rain0*ones(1,Ny+1); %moor
Rr(1:Ny+1) = nrain*Rain0*ones(1,Ny+1); %res
nrainc = 1;
ncc = 0;
%% Video?

if (VIDEO == 1)
    % make directory path
    outdirv = strcat(pwd, {'/'}, 'mov/');
    outdirv = strjoin(outdirv);
    
    res = num2str(Nk);
    Tmax = strrep(num2str(Tend), '.', '_');
    
    outnamev = sprintf('wetro_Nk=%s_Tend=%s',res,Tmax);
    outnamev = strcat(outdirv,outnamev);
    v = VideoWriter(outnamev);
    v.FrameRate = 10;
    v.Quality = 50;
    open(v);
    
end

%% Step forward in time...
CFL = 0.35; % CFL number for stable time-stepping
while (tijd  <= Tend) % All simple explicit time stepping
    % Randomised rainfall 
    % (3, 7, 5, 1)/16; [0,3]/16, [3,10]/16, [10,15]/16, [15,16]/16
    switch nrainc
        case 1
            while (tijd >= tunit)
                
                ncc = ncc+1;
                tunit  = tunit+TimeUnit;
                trand1 = rand;
                trand2 = rand;
                
                %%% "galton board 1": RAINFALL AMOUNT
                if (trand1 < 3/16)
                    nrain = 1; % Rain0
                else
                    if (trand1 < 10/16)
                        nrain = 2; % 2*Rain0
                    else
                        if (trand1 < 15/16)
                            nrain = 4; % 4*Rain0
                        else
                            nrain = 9; % 9*Rain0
                        end
                    end
                end
                
                %%% "galton board 2": RAINFALL LOCATION
                if (trand2 < 3/16)
                    nloc = 1; % Reservoir
                    Rm(1:Ny+1) = 0*Rain0*ones(1,Ny+1); %moor
                    Rr(1:Ny+1) = nrain*Rain0*ones(1,Ny+1);%res
                    not = nrain;
                else
                    if (trand2 < 10/16)
                        nloc = 2; % Moor and reservoir
                        Rm(1:Ny+1) = nrain*Rain0*ones(1,Ny+1); %moor
                        Rr(1:Ny+1) = nrain*Rain0*ones(1,Ny+1); %res
                        not = 2*nrain;
                    else
                        if (trand2 < 15/16)
                            nloc = 3; % Moor
                            Rm(1:Ny+1) = nrain*Rain0*ones(1,Ny+1); %moor
                            Rr(1:Ny+1) = 0*Rain0*ones(1,Ny+1);%res
                            not = nrain;
                        else
                            nloc = 4; % No rain in catchment
                            Rm(1:Ny+1) = 0*Rain0*ones(1,Ny+1); %moor
                            Rr(1:Ny+1) = 0*Rain0*ones(1,Ny+1); %res
                            not = 0;
                        end
                    end
                end
                
%                 switch not
%                     case 0
%                         nc1 = nc1+1;
%                     case 2
%                         nc2 = nc2+1;
%                     case 4
%                         nc4 = nc4+1;
%                     case 8
%                         nc8 = nc8+1;
%                     case 16
%                         nc16 = nc16+1;
%                 end
                
                nchisto(ncc) = not;
                fprintf('tunit: %g nrain: %d noloc: %d not:%d\n',tunit,nrain,nloc,not);
                
            end
            
        case 2
            
            if (tijd > Train)
                Rm(1:Ny+1) = Rain0*ones(1,Ny+1);  % 1mm/s
                Rr(1:Ny+1) = Rain0*ones(1,Ny+1); % 1mm/s
            else
                Rm(1:Ny+1) = Rain0*ones(1,Ny+1);
                Rr(1:Ny+1) = Rain0*ones(1,Ny+1); % 1mm/s
            end
            
    end
    %%
    
    % River: determine hydraulic radius, h and dh/dA for time step restrictions
    % and numerical flux 
    [h(1), dhdA(1)] = xsec_hAs(U(1,1),0.5*(-Kk+0),geom,chan);
    [h(Nk+2), dhdA(Nk+2)] = xsec_hAs(U(1,Nk+2),0.5*(L+L+Kk),geom,chan);
    [area(1), Wp(1), Rh(1)] = xsec_Ahs(h(1),0.5*(-Kk+0),geom,chan);
    [area(Nk+2), Wp(Nk+2), Rh(Nk+2)] = xsec_Ahs(h(Nk+2),0.5*(L+L+Kk),geom,chan); 
        
    for j = 2:Nk+1
        [h(j), dhdA(j)] = xsec_hAs(U(1,j),0.5*(s(j-1)+s(j)),geom,chan);
        [area(j), Wp(j), Rh(j)] = xsec_Ahs(h(j),0.5*(s(j-1)+s(j)),geom,chan);
    end
    
    %moor
    hmg = max(0.5*(hm(3:Ny+1)+hm(2:Ny)+0.5*(hm(2:Ny)+hm(1:Ny-1))));
    hmg = max(hmg,hm(Ny)+hm(Ny+1));
    
    %kinematic
    Vr = (wr*hr(1:Nx)./(2*hr(1:Nx)+wr)).^(2/3)*sqrt(-dbdxmean)/Cm; % kinematic velocity
    
    dtr = dx/max(abs(Vr)); %dt kinematic river
    dtm = dy^2/(alph*g*max(hmg,0.001)); %dt moor?
    dtAu = min(Kk./max(abs(U(2,:)./U(1,:) - sqrt(geom.g*U(1,:).*dhdA)),...
        abs(U(2,:)./U(1,:) + sqrt(geom.g*U(1,:).*dhdA)))); %dt st venant
    dt = CFL*min(dtr,min(dtm,dtAu)); % take min timestep with CFL
    
    tijd = tijd+dt;
    
    if (tijd > tmeas)
        dt = dt - (tijd - tmeas) + 1.e-10;
        tijd = tmeas+1.e-10;
    end
    
    nt=nt+1;
    %
    %% Update canal sections
    %
    h3co = h3c; 
    h2co = h2c; 
    h1co = h1c;
    tesfac = 0.5*(1+sign(h1co-hr(nxLc1)));
    tesfac = 1;
    % Q3c, Q2c, Q1c
    Qc3 = wc*sqrt(g)*Cf*max(h3co-Pw3,0)^(3/2);        % units: m m^(1/2)/s m^(3/2) = m^3/s
    Qc2 = wc*sqrt(g)*Cf*max(h2co-Pw2,0)^(3/2);        % units: m m^(1/2)/s m^(3/2) = m^3/s
    Qc1 = wc*sqrt(g)*Cf*max(h1co-Pw1,0)^(3/2)*tesfac; % units: m m^(1/2)/s m^(3/2) = m^3/s
    % first order approximations
    h3c = h3co+(dt/(Lsec3*wc))*( gam*mpor*sigme*wv*alph*g*0.5*(hm(2)^2-h3co^2)/dy-Qc3 ); % corrected  Jan 2017: needs correction h2co -> h3co
    h2c = h2co+(dt/(Lsec2*wc))*( Qc3-Qc2 );
    h1c = h1co+(dt/(Lsec1*wc))*( Qc2-Qc1 );
    tes = 0.5*(1+sign(bx(1:Nx)+hr(1:Nx)-Hcxb(1:Nx)-h1co));
    % h1c = h1c+nswitch*(dt/(Lc1*wc))*Cf*sqrt(g)*max(bx(1:Nx)+hr(1:Nx)-Hcc1,0).^(3/2).*(1+sign(xxm-Lc2)).*tes*0.5*transpose(Kk); % addition of river water wrong
    %
    %% Update groundwater model (moor)
    %
    num = dt/dy^2;
    %
    % first grid point
    % Option a) canal bottom coincide with Hele-Shaw moor bottom such that hm(1) = hc1; b) hm(1) = hm(1)+num;
    hmo = hm;
    hm(1) = h3co; % Jan 2017: this is wrong needs correction h2co -> h3co
    % interior points
    hm(2:Ny) = hmo(2:Ny)+num*alph*g*(0.5*(hmo(3:Ny+1)+hmo(2:Ny)).*(hmo(3:Ny+1)-hmo(2:Ny))-0.5*(hmo(2:Ny)+hmo(1:Ny-1)).*(hmo(2:Ny)-hmo(1:Ny-1)))+dt*Rm(2:Ny)/(mpor*sigme);
    % wall at last grid point so no flux
    hm(Ny+1) = hmo(Ny+1)+num*alph*g*(hmo(Ny)+hmo(Ny+1))*(hmo(Ny)-hmo(Ny+1))+dt*Rm(Ny+1)/(mpor*sigme);
    
    %
    %% Update reservoir hres(t)
    %
    hreso = hres;
    Qresw = Cf*wres*sqrt(g)*max(hreso-Pwr,0)^(3/2);
    hres = hreso+dt*Rr(1)-(dt/(Lyy*wres))*Qresw;
    %
     %% Update river kinematic
%     %
%     Ar0 = hr0*wr; % in m hrofunc(tijd,wr);
%     Aro = Ar;
%     Qh12(1) = Ar0*(hr0*wr/(2*hr0+wr))^(2/3)*sqrt(-dbdxmean)/Cm;
%     Qh12(2:Nx+1) = Aro(1:Nx).*(wr*hr(1:Nx)./(2*hr(1:Nx)+wr)).^(2/3)*sqrt(-dbdxmean)/Cm; % upwind
%     Ar(1:Nx) = Aro(1:Nx)-dt*(Qh12(2:Nx+1)-Qh12(1:Nx))./Kk;
%     % Flooding river terms
%     Ar(nxLv) = Ar(nxLv) + dt*((1-gam)*mpor*sigme*wv*alph*g*0.5*(hmo(2)^2-h3co^2)/dy)/Kk(nxLv); % addition of moor water corrected  Jan 2017: needs correction h2co -> h3co
%     Ar(nrer) = Ar(nrer) + dt*Qresw/Kk(nrer);                                                   % addition of reservoir water
%     Ar(nxLc1)= Ar(nxLc1)+ nswitch2*dt*Qc1/Kk(nxLc1);                                           % addition of canal water
%     %
%     % Ar(1:Nx) = Ar(1:Nx) - nswitch*dt*Cf*sqrt(g)*max(bx(1:Nx)+hr(1:Nx)-Hcc1,0).^(3/2).*(1+sign(xxm-Lc2))*0.5.*tes; % addition of overflown canal bank wrong
%     %
%     hr = Ar/wr;
%     Vr = (wr*hr(1:Nx)./(2*hr(1:Nx)+wr)).^(2/3)*sqrt(-dbdxmean)/Cm; % velocity
%     Qr = Vr.*Ar; % dscharge
    %
    %
    %% update river Au 
    for j = 1:Nk+1
        [Flux(:,j), SL(j), SR(j), VNC(:,j)] = NCPflux_Au(U(:,j),U(:,j+1),s(j),geom,chan);
    end
    
    %%%----- Determine hydraulic radius, h and dh/dA -----%%%
    [h(1), dhdA(1)] = xsec_hAs(U(1,1),0.5*(-Kk+0),geom,chan);
    [h(Nk+2), dhdA(Nk+2)] = xsec_hAs(U(1,Nk+2),0.5*(L+L+Kk),geom,chan);
    [area(1), Wp(1), Rh(1)] = xsec_Ahs(h(1),0.5*(-Kk+0),geom,chan);
    [area(Nk+2), Wp(Nk+2), Rh(Nk+2)] = xsec_Ahs(h(Nk+2),0.5*(L+L+Kk),geom,chan); 
        
    for j = 2:Nk+1
        [h(j), dhdA(j)] = xsec_hAs(U(1,j),0.5*(s(j-1)+s(j)),geom,chan);
        [area(j), Wp(j), Rh(j)] = xsec_Ahs(h(j),0.5*(s(j-1)+s(j)),geom,chan);
    end
    
%     %%%----- compute topographic terms as per Audusse -----%%%
%     
%     Sb(2,:) = 0.5*g*(Uminus(1,2:end).^2 - Uplus(1,1:end-1).^2);
    
    %%%----- compute extraneous forcing terms S(U) -----%%%
%     Rh = 0.009;
    
    S(1,:) = 0;
    S(2,:) = -geom.g*U(1,:)*geom.dbds - geom.g*geom.Cm^2*U(2,:).*abs(U(2,:)./U(1,:))./Rh.^(4/3);

%     %%%---- determine timestep for stability using wave eigen-speeds -----%%% 
%     
%     dt = cfl*min(Kk./max(abs(U(2,:)./U(1,:) - sqrt(geom.g*U(1,:).*dhdA)),...
%         abs(U(2,:)./U(1,:) + sqrt(geom.g*U(1,:).*dhdA))));
%     
%     %%%----- update time given new time step -----%%%
%     tn = tn + dt;
%     
%     if (tn > tmeasure)
%         dt = dt - (tn - tmeasure) + 1.e-10;
%         tn = tmeasure+1.e-10;
%     end
    
    %%% P fluxes as per the NCP theory
    Pp = 0.5*VNC + Flux;
    Pm = -0.5*VNC + Flux;
    
    %%%----- integrate forward to next time level -----%%% 
    
    if (BC == 1) % periodic
        
%         UU = U - dt*(Pp(:,2:Nk+1) - Pm(:,1:Nk))./Kk + dt*Sb./Kk + dt*S;
        UU = U - dt*(Pp(:,2:Nk+1) - Pm(:,1:Nk))./Kk + dt*S; % NOTE: not updated

    elseif (BC == 2) % neumann
        
        %interior
        UU(:,2:end-1) = U(:,2:end-1) - dt*(Pp(:,3:Nk) - Pm(:,2:Nk-1))./Kk + dt*S(:,2:end-1);
        %ghosts
        UU(:,1) = UU(:,2);
        UU(:,Nk) = UU(:,Nk-1);
        
    elseif (BC == 3) % specified inflow
       
        % interior
        UU(:,2:end-1) = U(:,2:end-1) - dt*(Pp(:,2:end) - Pm(:,1:end-1))./Kk + dt*S(:,2:end-1);
        UU(:,2:end-1) = U(:,2:end-1) - dt*(Pp(:,2:end) - Pm(:,1:end-1))./Kk + dt*S(:,2:end-1);

        % ghosts
        UU(1,1) = U(1,2);
%         UU(2,1) = U0(2,1)+0.00005*sin(tn/(4*pi)); % Au: sine wave
%         UU(2,1) = U0(2,1) + 0.0004*exp(-((tijd-0.25*Tend).^2)/50); % Au: exp pulse
        UU(2,1) = U0(2,1);
        UU(:,end) = UU(:,end-1);

    end
    
    % Area: update point source terms S_A
    UU(1,nxLv) = UU(1,nxLv) + dt*((1-gam)*mpor*sigme*wv*alph*g*0.5*(hmo(2)^2-h3co^2)/dy)/Kk; % moor inflow
    UU(1,nrer) = UU(1,nrer) + dt*Qresw/Kk; % res inflow 
    UU(1,nxLc1)= UU(1,nxLc1)+ nswitch2*dt*Qc1/Kk; % canal 1 inflow
    
    % Discharge: update point source term --- u*S_A
    UU(2,nxLv) = UU(2,nxLv) + (UU(2,nxLv)/UU(1,nxLv))*dt*((1-gam)*mpor*sigme*wv*alph*g*0.5*(hmo(2)^2-h3co^2)/dy)/Kk; % moor inflow 
    UU(2,nrer) = UU(2,nrer) + (UU(2,nrer)/UU(1,nrer))*dt*Qresw/Kk; % res inflow
    UU(2,nxLc1)= UU(2,nxLc1)+ (UU(2,nxLc1)/UU(1,nxLc1))*dt*nswitch2*Qc1/Kk; % canal 1 inflow
%     
    %%%----- update arrays for A, Au and h -----%%%
    U = UU;    
    %ghosts
    h(1) = xsec_hAs(U(1,1),0.5*(-Kk+0),geom,chan);
    h(end) = xsec_hAs(U(1,end),0.5*(L+L+Kk),geom,chan);
    %interior
    for j = 2:Nk+1
        h(j) = xsec_hAs(U(1,j),0.5*(s(j-1)+s(j)),geom,chan);
    end
    
    %% Live plotting at intervals:
    %
    % if (mod(nt,nmeas)==0)
    if (tijd >= tmeas)
        
        tmeas = tmeas+dtmeas;
        % nt
        
        A = U(1,:);
        Au = U(2,:);
              
        fp = index_fp(50);
        ct = index_city(5);
        
        % NOTE: legends are SLOW for live plotting
        fs = 14; %fontsize
        mainfig = figure(111);
        mainfig.WindowState = 'maximized';
%         subplot(2,2,1);
%         subplot(2,6,2); 
        subplot(2,5,3);
        plot(tijd,h1c,'ok','linewidth',1); hold on;
        plot(tijd,h2c,'*c','linewidth',1); hold on;
        plot(tijd,h3c,'xb','linewidth',1); hold on;
        plot(tijd,hres/10,'*r','linewidth',1); hold on;
%         legend({'Canal 1','Canal 2', 'Canal 3', 'Res/10'},'Location','northwest','fontsize',fs);
        ylim([0,0.025]);
        xlabel('t [s]','fontsize',fs);
%         ylabel('h_{1c}(t) (black), h_{3c}(t) (blue), h_{3c}(t) (cyan) h_{resr}(t)/10','linewidth',1); hold on;
        ylabel('Water depth h [m]', 'fontsize',fs); hold on;
        
%         subplot(2,2,2);
%         subplot(2,6,4);
        subplot(2,5,8);
        plot([tijdo,tijd],[hro,h(1)],'b','linewidth',1); hold on; %s=0
        plot([tijdo,tijd],[hrLo,h(nxLc1+2)],'k','linewidth',1); hold on; %s=L_c1
        plot([tijdo,tijd],[0.02,0.02],':r','linewidth',1); hold on; %h flood?
%         axis([0 Tend 0 max(h2c)]);
%         legend({'h(s=0,t)','h(s=L_{1c},t)'},'Location','southeast','fontsize',fs);
        ylim([0,0.025]);
        xlabel('t [s]', 'fontsize',fs);
        ylabel('River level h [m]', 'fontsize',fs); hold on;
        
%         subplot(2,2,3);
%         subplot(2,6,3);
        subplot(2,5,2);        
        plot(yy,hm,'-','linewidth',2);
%         legend({'Groundwater level'},'Location','southeast','fontsize',fs);
        ylim([0,0.12]);
        xlabel('y [m]', 'fontsize',fs);
        ylabel('h_m(y,t) [m]', 'fontsize',fs);
        
%         subplot(2,2,4);
%         subplot(2,6,1);
        subplot(2,5,1);
        plot(tijd,Rm(1)/Rain0,'oc','linewidth',1); hold on;
        plot(tijd,Rr(1)/Rain0,'xb','linewidth',1); hold on;
        plot(tijd,(Rm(1)+Rr(1))/Rain0,'*r','linewidth',1); hold  on;
%         legend({'Moor','Res', 'Both'},'Location','northeast','fontsize',fs);
        ylim([0,18]);
        xlabel('t [s]', 'fontsize',fs);
        ylabel('Rainfall factor', 'fontsize',fs);
        %
%         sgtitle(['Wetropolis flood and rainfall demonstrator: t = ',num2str(tijd),' s'])

        %
%         figure(14);
%         %
%         subplot(2,1,1);
% %         plot(xxm,Vr/5,'-k','linewidth',2); hold on;
%         plot(xxm,h(2:end-1),'-b','linewidth',2); hold off;
%         xlabel('x (m)');
%         ylabel('h_r(x,t) (blue), V_r(x,t)/5 (black)');
%         ylim([0,0.025]);
%         %
%         subplot(2,1,2);
%         plot(xxm,bx,'-k','linewidth',2); hold on % river topog
%         plot(xxm,B(2:end-1),'-k','linewidth',2); hold on % river topog
%         plot(xxm(1:nxLc1),Hcx(1:nxLc1),'-.r','linewidth',1);  hold on % canal berm?
%         plot(xxm(1:nxLc1),Hcxb(1:nxLc1),'-.k','linewidth',2);  hold on %canal topog
%         Hcxb12 = (Hcc3+h3c)*0.5*(1+sign(Lc3-xxm(1:nxLc1)))+(Hcc2+h2c)*0.5*(1+sign(xxm(1:nxLc1)-Lc3)).*(1+sign(Lc2-xxm(1:nxLc1)))*0.5+(Hcc1+h1c)*0.5*(1+sign(xxm(1:nxLc1)-Lc2))-canalmaxdepth;
%         plot(xxm(1:nxLc1),Hcxb12(1:nxLc1),'-.b','linewidth',1);  hold on %canal depth?
%         plot(xxm,B(2:end-1)+h(2:end-1),'-b','linewidth',2); hold off % river level
%         xlabel('x [m]', 'fontsize',fs);
%         ylabel('Height z [m]', 'fontsize',fs);
%         ylim([0,0.07]);
        %
        
%         fg = figure(114);
%         subplot(2,2,1);
%         subplot(2,6,[11 12]);
        subplot(2,5,9);
        for k = 1:(length(s)-1)
            plot([s(k), s(k+1)],[h(k+1),h(k+1)],'b-'); hold on;
        end
        plot([chan.s_r, chan.s_r],[0,0.04],'r:'); hold on;
        plot([chan.s_m, chan.s_m],[0,0.04],'r:'); hold on;
        plot([Lc1, Lc1],[0,0.04],'r:'); hold on;
        plot([s(fp), s(fp)],[0,0.04],'k:'); hold on;
        plot([s(ct), s(ct)],[0,0.04],'k:'); hold on;
        fill([chan.LR1, chan.LR2,chan.LR2,chan.LR1],[0,0,geom.hc,geom.hc],...
            'r','FaceAlpha',0.1,'LineStyle','none');
        hold off;
%         text(0.85*xmax,0.9*hmax,['t=',num2str(tn)]);
        axis([0 L 0.01 0.03]);
        xlabel('s','fontsize',fs); ylabel('h(s,t)','fontsize',fs);
        title([]);

%         subplot(2,2,2);
%         subplot(2,6,[5 6]);
        subplot(2,5,10);
        for k = 1:(length(s)-1)
            plot([s(k), s(k+1)],[Au(k+1),Au(k+1)],'b-'); hold on;
        end
        plot([chan.s_r, chan.s_r],[0,0.04],'r:'); hold on;
        plot([chan.s_m, chan.s_m],[0,0.04],'r:'); hold on;
        plot([Lc1, Lc1],[0,0.04],'r:'); hold on;
        hold off;
        axis([0 L 0.0001 0.0004]);
        xlabel('s','fontsize',fs); ylabel('Au(s,t)','fontsize',fs);
        title([]);
        
%         subplot(2,2,4);
%         subplot(2,6,10);
        subplot(2,5,7);
        [X,Y,Xc,Yc,~] = plot_xsec_hAs(A(ct+1),s(ct),geom,chan);
        plot(Xc,Yc,'k', 'Linewidth',2); hold on;
        fill(X,Y,'b','FaceAlpha',0.1); hold off;
        text(Xc(1)+0.01,0.9*Yc(1),['t=',num2str(tijd)],'fontsize',fs,'HorizontalAlignment', 'left');
        text(Xc(1)+0.01,0.8*Yc(1),['s=',num2str(s(ct))],'fontsize',fs,'HorizontalAlignment', 'left');
        axis([min(Xc) max(Xc) min(Yc) max(Yc)]);
        
%         subplot(2,2,3);
%         subplot(2,6,9);
        subplot(2,5,6);
        [X,Y,Xc,Yc,hc] = plot_xsec_hAs(A(fp+1),s(fp),geom,chan);
        plot(Xc,Yc,'k', 'Linewidth',2); hold on;
        fill(X,Y,'b','FaceAlpha',0.1); hold off;
        text(Xc(1)+0.01,0.9*Yc(1),['t=',num2str(tijd)],'fontsize',fs,'HorizontalAlignment', 'left');
        text(Xc(1)+0.01,0.8*Yc(1),['s=',num2str(s(fp))],'fontsize',fs,'HorizontalAlignment', 'left');
        axis([min(Xc) max(Xc) min(Yc) max(Yc)]);
        
%         figure(13);
%         subplot(1,2,1);
%         subplot(2,6,7);
        subplot(2,5,4);
%         ntot = Tend/TimeUnit;
        histogram(nchisto, 'Normalization','pdf');
        xlabel('Rain / wd', 'fontsize',fs);
        ylabel('Density', 'fontsize',fs);
        axis([0 20 0 0.4]);
%         subplot(1,2,2);

%         subplot(2,6,8);
        subplot(2,5,5);
        plot(h(nxLc1+2), U(2,nxLc1+2),'xk'); hold on;
        xlabel('h', 'fontsize',fs);
        ylabel('Q', 'fontsize',fs);
        axis([0.01 0.03 0.0001 0.0003]);
        drawnow; pause(0.1);
        
        if (VIDEO == 1)
            frame = getframe(mainfig);
            writeVideo(v,frame);
        end

        hro = h(1);
        hrLo = h(nxLc1+2);
        tijdo = tijd;
        
    end
    
end

if (VIDEO == 1)
    close(v);
end
% figure(13);
% ntot = Tend/TimeUnit;
% histogram(nchisto);

