%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       WETROPOLIS RAINFALL AND FLOOD DEMONSTRATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - RIVER MODEL: ST. VENANT EQ.
% - GROUNDWATER MODEL: further developed (See below)
% - RESERVOIR MODEL: now coupled to canal and river 
% - CANALS: correct moor coupling (optional), res coupling (optional).

% v3 builds on run_wetro_2020vr.m. Updates here:

% > clean up code: lots of redundant parts and repetition. DONE/ongoing.
% > Advanced groundwater model (FEM) with connecting channel. PARTIALLY
% DONE, needs checking.
% > save rainfall data for repeat experiments 
% > BCs: steep bedslope at boundaries to enforce inflow/outflow via
% numerical flux?
% > EnKF for river model?
% > simple control problem: constrain river depth less than hT = 0.02 in
% city by, e.g., turning off reservoir outflow (i.e., raise the weir height so reservoir acts as storage
% buffer)

clear;
nparaflag = 2; % 1 is Onno's original design case; 2 is Luke/Onno/real case of parameter choice'
VIDEO = 0; % 0 to run without saving a vid, 1 to save.
RAINSAVE = 0; % 1 to save randomly-generated rainfall, 0 to not save
RAINLOAD = 0; % 1 to load saved rainfall (must specify file name), 0 
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

gam_m = 0.0; % proportion of water entering canal from moor (zero in physical model)
gam_r = 0.2; % proportion of water entering canal from reservoir (unknown) 

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
Cm = 0.02;
% Cm = 0.06; % Manning coefficient, see Table in Munson and or wiki
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
%% (b) St Venant river model
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

% test steep bed at s=0,L for effieicent boundaries.

%%%----- channel length parameters ----%%%
chan.LR3 = 5.2; % total length of channel
% chan.LR3 = Lx; % total length of channel

chan.LR1 = 3.8; % length upto city region
chan.LR2 = 4.2; % length from end of city region

chan.LR11 = 3.6;% transition zone from fp to c [LR11, LR1]
chan.LR22 = 4.4;% transition zone from c to fp [LR2, LR22]

chan.s_r = 0.932; % reservoir influx
chan.s_m = 2.038; % moor influx
% sr = chan.s_r;
% sm = chan.s_m;

chan.tr = 50; % severity of transition 

%
%%%----- # of equations in system -----%%%
Neq = 2; % U = (A, Au)

%%%----- Set up grid with curvilinear coord s -----%%%
L=chan.LR3; %length of domain
Nk=500; %number of gridcells (excluding ghost)
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

% test steep bed at s=0,L for inflow/outflow boundaries?
dbds_vec = geom.dbds*ones(1,Nk+2);
jj = 3; % how many cells to steepen at either end?
dbBC = -10.0; % steepness?
dbds_vec(1:jj) = dbBC; 
dbds_vec(end-jj+1:end) = dbBC;

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
        sm = 0.932; % m; % inflow of water off the Hele-Shaw moor in m
    case 2 % Luke/real settings
        hr0 = 0.0135; % m
        Ly = 0.925; % m
        wv = 0.095; % width Hele-Shaw moor cell in m
        sm = 2.038; % m; % inflow of water off the Hele-Shaw moor in m
end

%% (1) finite diff or (2) finite element model
moormodel = 1;
switch moormodel
    
    case 1 %FD
        
        Ny = 20;
        dy = Ly/Ny;
        yy = 0:dy:Ly;  % Regular FD or FV grid
        hm0 = 0*yy;     % Define moor depth variable on grid
        sigm = 0.8;    % fraction of moor pores filled
        sigm0 = 0.1;   % fraction of water remaining in moor pores after water has seeped out
        sigme = sigm;  %
        mpor = 0.3;    % porosity moor
        nu = 10^(-6);  % viscosity water m^2/s
        kperm = 10^-8; % permeability
        alph = kperm/(mpor*nu*sigm);
        
        hcm = 0.0;
        Pwm = 0.02;
        Lc = 0.1;
        
    case 2 %FEM
        
        Lc = -0.05; % connecting canal [m]
        Lca = abs(Lc);
        sigm = 0.8;    % fraction of moor pores filled
        sigm0 = 0.1;   % fraction of water remaining in moor pores after water has seeped out
        sigme = sigm;  %
        mpor = 0.3;    % porosity moor
        nu = 10^(-6);  % viscosity water m^2/s
        kperm = 10^-8; % permeability
        alph = kperm/(mpor*nu*sigm);
        %
        Ny = 20;   % number of elements
        Nn = Ny+1; % number of nodes
        Mm = zeros(Nn,Nn);
        Sm = zeros(Nn,Nn);
        
        dy = Ly/Ny;
        yy = transpose(0:dy:Ly);  % Regular FEM grid
        Kky = (yy(2:Nn)-yy(1:Nn-1)); % finite element lengths
        nini = 1;
        g3 = sqrt(1.0/3.0);
     
        
        hcm = 0.0;
        Pwm = 0.02;
        Lc = 0.1;
        
        bini = zeros(Nn,1);
        for nk=1:Ny
            for na = 1:2
                ii = (nk-1)+na;
                bini(ii) = bini(ii)+bgina(Kky(nk),yy(nk),yy(nk+1),na,g3,nini,Ly);
                for nb = 1:2
                    jj = (nk-1)+nb;
                    Mm(ii,jj) = Mm(ii,jj)+Mmalbe(Kky(nk),yy(nk),yy(nk+1),na,nb,g3);
                    Sm(ii,jj) = Sm(ii,jj)+Smalbe(Kky(nk),yy(nk),yy(nk+1),na,nb,g3);
                end
            end
        end
        
        hm0 = zeros(Nn,1);
        hm0(1:Nn) = Mm\bini;
        
        Ml = Mm;
        Ml(1,1) = Ml(1,1)+Lca/(mpor*sigme);
        
end

% moor canal = 0
%% Reservoir
% Reservoir Hele-Shaw cell with adjustable weir
%
switch nparaflag
    case 1
        sr = 2.038;  % m
        wres = 0.01;   % m	%
        Lyy  = 0.4;    % length reservoir in
    case 3
        sr = 2.038;  % m
        wres = 0.01;   % m	%
        Lyy  = 0.4;    % length reservoir in
    case 2
        sr = 0.932;  % m
        wres = 0.123;  % m	%
        Lyy  = 0.293;    % length reservoir in m
end
Pwr  = 0.10; % weir height


%% Moor, res, canal 1 coupling location
%
%
nxsr = ceil(sr/Kk)+1; %  grid box in which res water is added in the river
nxsm = ceil(sm/Kk)+1;    % grid box in which moor water is added in the river
nxLc1 = ceil(Lc1/Kk)+1;  % grid box for regular grid in which canal water segment one is added
%
% nxr = index_r(1)-1; %  grid box in which res water is added in the river
% nxsm = index_m(1)-1;    % grid box in which moor water is added in the river
% nxLc1 = ceil(Lc1/dx);  % grid box for regular grid in which canal water segment one is added

%% Rainfall and time
switch nparaflag
    case 1 % Onnos original
        Rain0 = 0.003; % 0.002 works m/s
        Rain = Rain0;
    case 2
        Rain0 = 1.5*0.00013656; % This is r0 Eq. 17 in HESS article
        Rain = Rain0;
end
%
TimeUnit = 10; % wd % 9 hrs from Kildwick till Armley times 3 is 27hrs; 12-20s across
Train = Lc2*wc*0.01/(Ly*wv*Rain); % TK: what is Train? unused (in case 2 below)
%
%
Vrate = Ly*wv*Rain0*[1,2,4,8,18]*TimeUnit*1000; % flow rates in liters/TimeUnit (Eq. 18 in HESS?)
%
rainfac = [0,1,2,4,8,9,18];
rainpdf = [16,24,77,89,35,8,7]/256;

tplot = 100;
Tend = tplot*TimeUnit;
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

ini = 2;
switch ini
    case 1 % Mass check starting from absolute drought conditions and then it starts raining
        % canal sections
        h1c = 0;
        h2c = 0;
        h3c = 0;
        % River
        Ar = hr0*wr*ones(1,Nx);
        U = U0;
        h = h0;
        % Hele-Shaw moor
        hm = hm0;
        % reservoir
        hres = 0; 
    case 2 % canals full to weirs, reservoir empty, groundwater level with river level at sm
        % canal sections
        h1c = Pw1;
        h2c = Pw2;
        h3c = Pw3;
        % River
        Ar = hr0*wr*ones(1,Nx);
        U = U0;
        h = h0;
        % Hele-Shaw moor
%         hm = hr0*ones(1,Ny+1); % initial river level
        hm = hm0;
        % reservoir
        hres = 0; 
end

hro = h(1);
hrLo = h(nxLc1);
tijdo = tijd;
%

nt = 0;
nswitch = 0;
nswitch2 = 1;
tunit = 0; % why tunit = 1 originally?
nrain = 4; 

% Rm(1:Ny+1) = 0*nrain*Rain0*ones(1,Ny+1); % why Rm and Rr non-zero
% initially?
% Rr(1:Ny+1) = 0*nrain*Rain0*ones(1,Ny+1); %

% nchisto = zeros(1,tplot);
rainm_save = zeros(Ny+1,tplot);
rainr_save = zeros(Ny+1,tplot);

nrainc = 1;
ncc = 0;
%% Video? Save?

if (RAINSAVE == 1)
    
    % make directory path
    outdirs = strcat(pwd, {'/'}, 'raindata/');
    outdirs = strjoin(outdirs);
    
    Tstr = strrep(num2str(Tend), '.', '_');
    
    outnames = sprintf('rain#1_tmax=%s',Tstr);
    
end

if (RAINLOAD == 1)
    
    rain = load('raindata/rain#1_tmax=1000.mat');
    Rmload = rain.raindata.rainm;
    Rrload = rain.raindata.rainr;    
    
    nrainc = 2;
end

if (VIDEO == 1)
    % make directory path
    outdirv = strcat(pwd, {'/'}, 'mov/');
    outdirv = strjoin(outdirv);
    
    res = num2str(Nk);
    Tstr = strrep(num2str(Tend), '.', '_');
    
    outnamev = sprintf('wetro5_Nk=%s_Tend=%s',res,Trt);
    outnamev = strcat(outdirv,outnamev);
    v = VideoWriter(outnamev);
    v.FrameRate = 5;
    v.Quality = 50;
    open(v);
    
end

%% Step forward in time...
CFL = 0.3; % CFL number for stable time-stepping
while (tijd  <= Tend) % All simple explicit time stepping
    %% Rainfall
    switch nrainc
        case 1 %Randomised rainfall 
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
                
                rainm_save(:,ncc) = Rm/Rain0;
                rainr_save(:,ncc) = Rr/Rain0;
                
                nchisto(ncc) = not;
                fprintf('tunit: %g nrain: %d noloc: %d not:%d\n',tunit,nrain,nloc,not);
                
            end
        
        case 2 % using saved rain 
            
            while (tijd >= tunit)
                
                ncc = ncc+1;
                tunit  = tunit+TimeUnit;
                
                Rm = Rmload(:,ncc)/Rain0;
                Rr = Rrload(:,ncc)/Rain0;
                
            end
            
        case 3 % uniform light rainfall continuous in time
            
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
    dtm = 0.5*dy^2/(alph*g*max(hmg,0.001)); %dt moor
    dtAu = min(Kk./max(abs(U(2,:)./U(1,:) - sqrt(geom.g*U(1,:).*dhdA)),...
        abs(U(2,:)./U(1,:) + sqrt(geom.g*U(1,:).*dhdA)))); %dt st venant
    dt = CFL*min(dtr,min(dtm,dtAu)); % take min timestep with CFL
    
    tijd = tijd+dt;
    
    if (tijd > tmeas)
        dt = dt - (tijd - tmeas) + 1.e-10;
        tijd = tmeas+1.e-10;
    end
    
%     nt=nt+1;
    %
    
    %% Update reservoir: hres(t)
    %
    hreso = hres;
    Qresw = Cf*wres*sqrt(g)*max(hreso-Pwr,0)^(3/2);
    hres = hreso+dt*Rr(1)-(dt/(Lyy*wres))*Qresw;
    
    %% Update groundwater model (moor)
    %
    switch moormodel
        
        case 1 % FD
            num = dt/dy^2;
            %
            hmo = hm;
            
            % first grid point
%             hm(1) = h3co; %Canal-3 (?) Jan 2017: this is wrong needs correction h2co -> h3co -- why???
%             hm(1) = h(nxsm); %River at sm
            hm(1) = hcm;% Moor canal 

            % interior points
            hm(2:Ny) = hmo(2:Ny)+num*alph*g*(0.5*(hmo(3:Ny+1)+hmo(2:Ny)).*(hmo(3:Ny+1)-hmo(2:Ny))-0.5*(hmo(2:Ny)+hmo(1:Ny-1)).*(hmo(2:Ny)-hmo(1:Ny-1)))+dt*Rm(2:Ny)/(mpor*sigme);
            % wall at last grid point so no flux
            hm(Ny+1) = hmo(Ny+1)+num*alph*g*(hmo(Ny)+hmo(Ny+1))*(hmo(Ny)-hmo(Ny+1))+dt*Rm(Ny+1)/(mpor*sigme);
            % outflow
            Qmoor = 0.5*mpor*sigme*wv*alph*g*(hmo(2)^2-hm(1)^2)/dy;
            
            % connecting canal?
            hcmo = hcm;
            Qcm = wr*sqrt(g)*Cf*max(hcmo-Pwm,0)^(3/2);
            hcm = hcmo+(dt/(Lc*wv))*( Qmoor - Qcm );
            
            
        case 2 %FEM -- coupled but not behaving as expected
            hmo = hm;
            
            bini = zeros(Nn,1); % size(bini)
            
            for nk=1:Ny
                for na = 1:2
                    ii = (nk-1)+na;
                    bini(ii) = bini(ii)+rainterm(Kky(nk),yy(nk),yy(nk+1),na,g3,nini,Ly,1.0/(mpor*sigme),Rm(nk),alph*g,hmo(nk),hmo(nk+1));
                end
            end
            
            bini(1) = bini(1)-sqrt(g)*max(2*hcm/3,0.0)^(1.5)/(mpor*sigme);
            hm = Ml\(Ml*hmo+dt*bini);
            hcm = hm(1);
            
%             hmo(1) = hcm;
%             hmo(2:Nn) = transpose(hm(2:Nn));
            
            Qmoor = 0.5*mpor*sigme*wv*alph*g*(hm(2)^2-hm(1)^2)/dy;
            
            % connecting canal?
            hcmo = hcm;
            Qcm = wr*sqrt(g)*Cf*max(hcmo-Pwm,0)^(3/2);
            hcm = hcmo+(dt/(Lc*wv))*( Qmoor - Qcm );
    end
    %
    %
    %% Update canal sections: h3c(t), h2c(t), h1c(t)
    %
    h3co = h3c; 
    h2co = h2c; 
    h1co = h1c;
    
    % Q3c, Q2c, Q1c
    Qc3 = wc*sqrt(g)*Cf*max(h3co-Pw3,0)^(3/2);        
    Qc2 = wc*sqrt(g)*Cf*max(h2co-Pw2,0)^(3/2);        
    Qc1 = wc*sqrt(g)*Cf*max(h1co-Pw1,0)^(3/2); 
    
    % first order approximations
    h3c = h3co+(dt/(Lsec3*wc))*( gam_r*Qresw - Qc3 );
    h2c = h2co+(dt/(Lsec2*wc))*( Qc3 + gam_m*Qmoor - Qc2 ); % note: gam_m = 0 in real set-up
    h1c = h1co+(dt/(Lsec1*wc))*( Qc2-Qc1 );

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
    %% update river St Venant 
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
    
    %%%----- compute extraneous forcing terms S(U) -----%%%
    S(1,:) = 0;
    S(2,:) = -geom.g*U(1,:)*geom.dbds - geom.g*geom.Cm^2*U(2,:).*abs(U(2,:)./U(1,:))./Rh.^(4/3);
    S(2,:) = -geom.g*U(1,:).*dbds_vec - geom.g*geom.Cm^2*U(2,:).*abs(U(2,:)./U(1,:))./Rh.^(4/3);
    
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
        UU(1,1) = U0(1,1);
%         UU(2,1) = U0(2,1)+0.00005*sin(tn/(4*pi)); % Au: sine wave
%         UU(2,1) = U0(2,1) + 0.0004*exp(-((tijd-0.25*Tend).^2)/50); % Au: exp pulse
        UU(2,1) = U0(2,1);
        UU(:,end) = UU(:,end-1);

    end
    
    % Area: update point source terms S_A
    UU(1,nxsm) = UU(1,nxsm) + dt*(1-gam_m)*Qcm/Kk; % moor inflow
    UU(1,nxsr) = UU(1,nxsr) + dt*(1-gam_r)*Qresw/Kk; % res inflow 
    UU(1,nxLc1)= UU(1,nxLc1)+ nswitch2*dt*Qc1/Kk; % canal 1 inflow
    
    % Discharge: update point source term --- u*S_A
    UU(2,nxsm) = UU(2,nxsm) + (UU(2,nxsm)/UU(1,nxsm))*dt*(1-gam_m)*Qcm/Kk; % moor inflow 
    UU(2,nxsr) = UU(2,nxsr) + (UU(2,nxsr)/UU(1,nxsr))*dt*(1-gam_r)*Qresw/Kk; % res inflow
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
        
        mainfig = figure(111);
        mainfig.WindowState = 'maximized';
        live_plotting_panel;
        
        if (VIDEO == 1)
            frame = getframe(mainfig);
            writeVideo(v,frame);
        end

        hro = h(1);
        hrLo = h(nxLc1+2);
        tijdo = tijd;
        
    end
    
end

%% close vid and save rain?
if (VIDEO == 1)
    close(v);
end

if (SAVE == 1)
    
    raindata.rainm = rainm_save;
    raindata.rainr = rainr_save;
    
    save(fullfile(outdirs, outnames),'raindata');
    
end
