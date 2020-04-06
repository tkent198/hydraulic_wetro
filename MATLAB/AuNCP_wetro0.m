%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wetropolis Au model -- first attempts with 2 cross sections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FV/DG0 with NCP theory. Bed slope prescribed: db/ds = - const.

% 14/8/19

clear; 
%%
VIDEO = 0; % 1 = TRUE for save a video
SAVE = 0; % 1 = TRUE for save data

%%
%%%----- geometry  parameters ----%%%
geom.hr = 0.015; 
geom.wr = 0.05; % m width main river channel
geom.hf = 0.02; % m
geom.hf = 0.005;
geom.hc = 0.02;
geom.hc = geom.hr+geom.hf;
geom.wf = 0.1;  % m width flood plain
geom.wc = 0.1; 
geom.tana = geom.hf/geom.wf;

geom.g = 9.81;     % acceleration of gravity
geom.Cm = 0.02;    % Manning coefficient
geom.dbds = -0.01; % mean slope river bed


%% From WetropA1d.m: 

% Spatial set-up of Wetropolis channel: length LR3 = circa 5m
% chan.L1 = 0.780;
% chan.r2 = 0.125;
% chan.L2 = pi*r2;
% chan.L3 = 0.637;
% chan.r4 = 0.175;
% chan.L4 = pi*r4;
% chan.L5 = 0.637;
% chan.r6 = 0.125;
% chan.L6 = pi*r6;
% chan.L7 = 0.427;
% chan.L10 = 0.094; % +del8;
% chan.L8 = 0.027;
% chan.del8 = L10-L8;
% chan.L7 = L7-del8;
% chan.L8 = L8+del8;
% chan.L9 = 0.463; % -del8;
% chan.r11 = 0.175;
% chan.L11 = 0.5*pi*r11;
% chan.L12 = 0.420;
% chan.LR1 = L1+L2+L3+L4+L5+L6+L7;
% chan.LR2 = LR1+L8+L9;
% chan.LR3 = LR2+L10+L11+L12;

chan.LR3 = 4.2; % total length of channel

chan.LR1 = 3.4; % length upto city region
chan.LR2 = 3.8; % length from end of city region

chan.LR11 = 3.4;% transition zone from fp to c [LR11, LR1]
chan.LR22 = 3.8;% transition zone from c to fp [LR2, LR22]

chan.s_r = 0.932; %reservoir influx
chan.s_m = 2.038; %moor influx

chan.tr = 50; % severity of transition 

%%
%%%----- # of equations in system -----%%%
Neq = 2; % U = (A, Au)

%%%----- Set up grid -----%%%
L=chan.LR3; %length of domain
Nk=25*L; %number of gridcells (excluding ghost)
Nf=Nk+1; %number of nodes
Kk=L/Nk; %length of cell
s=0:Kk:L;%node location
sBC = [-Kk s L+Kk]; %node loc with ghosts

% locating floodplain/city
index_fp = [find(s < chan.LR1) find(s > chan.LR2)];
index_city = intersect(find(s > chan.LR1), find(s < chan.LR2));


cfl = 0.5; % CFL number for stable time-stepping

% BOUNDARY CONDITIONS
% Periodic BC = 1
% Neumann BC = 2
% specified inflow BC = 3
BC = 3; 

%% ----- Apply initial conditions ----- %%
[U0, B, h0] = initial_cond_wetro(s,Neq,Nk,geom,chan);

U0 = [U0(:,1) U0 U0(:,end)];
B = [B(1) B B(end)];
h0 = [h0(1) h0 h0(end)];

Z0 = B + h0;

% % check ICs: plot
% figure(1);
% for k = 1:(length(s)-1)
%     plot([s(k), s(k+1)],[B(k+1),B(k+1)],'b-'); hold on;
%     plot([s(k), s(k+1)],[Z0(k+1),Z0(k+1)],'k-'); hold on;
% end
% hold off;
% axis([0 L 0 0.06]);
% xlabel('s','fontsize',18); ylabel('h+b','fontsize',18);
% set(gcf, 'Position',  [100, 100, 1000, 600])


%% ----- time parameters -----%%
tmax = 100;
tn = 0;
Nmeas = 100;
dtmeasure = tmax/Nmeas;
tmeasure = dtmeasure;

%% TESTS: plots of water depth h in cross section A at location s: h = h(A(s,t),s)
% ind = index_fp(10); % index for position
% sx = s(ind);
% Ax = U0(1,ind);
% 
% Xcplot = zeros(index_city(1)-1,6);
% Ycplot = Xcplot;
% 
% for jj = 1:index_city(1)-1
%     [~,~,Xcplot(jj,:),Ycplot(jj,:),~] = plot_xsec_hAs(0,s(jj),geom,chan);
% end


% 
% while (tn < tmax)
%     tn = tn + dtmeasure;
%     [X,Y,Xc,Yc,~] = plot_xsec_hAs(Ax*(1+0.75*sin(tn/(8*pi))),sx,geom,chan);
%     figure(111);
%     plot(Xc,Yc,'k', 'Linewidth',2); hold on;
%     fill(X,Y,'b','FaceAlpha',0.1); hold off;
%     text(Xc(end),0.25*geom.hr,['Location, s = ',num2str(sx)],'fontsize',14,'HorizontalAlignment', 'right');
%     text(Xc(end),0.5*geom.hr,['Time, t = ',num2str(tn)],'fontsize',14,'HorizontalAlignment', 'right');
%     drawnow; pause(0.1);
% end




%% Define system arrays with ghost cells for BCs
Flux = zeros(Neq,Nk+1);
S = zeros(Neq,Nk+2); % for source terms (friction, bed slope, and rain)
UU = zeros(Neq,Nk+2);
SL = zeros(1,Nk+1); SR = zeros(1,Nk+1); % numerical speeds
VNC = zeros(Neq,Nk+1); % NCP integrals
area = zeros(1,Nk+2);
Wp = zeros(1,Nk+2);
Rh = zeros(1,Nk+2);
dhdA = zeros(1,Nk+2);

U = U0;
h = h0;

% indexing for reconstructing states as per Audusse
% left = [1:Nk+1]; 
% right = [2:Nk+2];

% set up arrays for storing data for hovmoller plots
hovz = zeros(Nmeas+1,Nk);
hovh = zeros(Nmeas+1,Nk);
hovA = zeros(Nmeas+1,Nk);
hovAu = zeros(Nmeas+1,Nk);

hovz(1,:) = Z0(2:end-1);
hovh(1,:) = h0(2:end-1);
hovA(1,:) = U0(1,(2:end-1)); 
hovAu(1,:) = U0(2,(2:end-1));


%% Video + save
if (SAVE == 1)
    
    % make directory path
    outdirs = strcat(pwd, {'/'}, 'data_AuNCP_wetro0/');
    outdirs = strjoin(outdirs);
    
    res = num2str(Nk);
    Tmax = strrep(num2str(tmax), '.', '_');
    
    outnames = sprintf('dat_Nk=%s_tmax=%s',res,Tmax);
    
end

if (VIDEO == 1)
    % make directory path
    outdirv = strcat(pwd, {'/'}, 'data_AuNCP_wetro0/');
    outdirv = strjoin(outdirv);
    
    res = num2str(Nk);
    Tmax = strrep(num2str(tmax), '.', '_');
    
    outnamev = sprintf('vid_Nk=%s_tmax=%s',res,Tmax);
    outnamev = strcat(outdirv,outnamev);
    v = VideoWriter(outnamev);
    v.FrameRate = 10;
    v.Quality = 100;
    open(v);
    
end

%% ---- Set up the numerical integration and plotting routines -----%%%
ihov = 2; % index for hov arrays
while (tn < tmax)
    
%     % State reconstructions for Audusse theory
%     h = U(1,:); hu = U(2,:);
%     hu(h < 1e-6) = 0; U(2,:) = hu;
%     
%     Bstar = max(B(left),B(right)); 
%     hminus = max(U(1,left) + B(left) - Bstar,0);
%     hplus = max(U(1,right) + B(right) - Bstar,0);
%     
%     ul = U(2,left).*spfun(@(x) 1./x, U(1,left)); 
%     ul = full(ul);
%     ul(isnan(ul)) = 0;% catch division by zero in dry regions and expand
%     uminus = ul;
%     ur = U(2,right).*spfun(@(x) 1./x, U(1,right)); 
%     ur = full(ur); % catch division by zero in dry regions and expand
%     ur(isnan(ur)) = 0;
%     uplus = ur;
%     
%     huminus = hminus.*uminus; huminus = [huminus(end), huminus];
%     huplus = hplus.*uplus; huplus = [huplus(end), huplus];
%     
%     hminus = [hminus(end), hminus];
%     hplus = [hplus(end), hplus];
%     uplus = [uplus(end), uplus];
%     uminus = [uminus(end), uminus];
%     
%     Bminus = [B(end), B(left)];
%     Bplus = [B(1), B(right)];
%     
%     
%     % reconstructed states
%     Uminus = [hminus; huminus; hrminus; hvplus];
%     Uplus = [hplus; huplus; hrplus; hvplus];
%    
    
    %%%----- Determine fluxes using flux function -----%%%
    
    % using ghost cells  
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

    %%%---- determine timestep for stability using wave eigen-speeds -----%%% 
    
    dt = cfl*min(Kk./max(abs(U(2,:)./U(1,:) - sqrt(geom.g*U(1,:).*dhdA)),...
        abs(U(2,:)./U(1,:) + sqrt(geom.g*U(1,:).*dhdA))));
    
    %%%----- update time given new time step -----%%%
    tn = tn + dt;
    
    if (tn > tmeasure)
        dt = dt - (tn - tmeasure) + 1.e-10;
        tn = tmeasure+1.e-10;
    end
    
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
        % ghosts
        UU(1,1) = U(1,2); % A -- is this OK??
%         UU(2,1) = U0(2,1)+0.00005*sin(tn/(4*pi)); % Au: sine wave
        UU(2,1) = U0(2,1) + 0.0004*exp(-((tn-0.25*tmax).^2)/50); % Au: exp pulse
        UU(:,end) = UU(:,end-1);

    end
    
    %%%----- update arrays for A, Au and h -----%%%
    U = UU;
    
    %ghosts
    h(1) = xsec_hAs(U(1,1),0.5*(-Kk+0),geom,chan);
    h(end) = xsec_hAs(U(1,end),0.5*(L+L+Kk),geom,chan);
    %interior
    for j = 2:Nk+1
        h(j) = xsec_hAs(U(1,j),0.5*(s(j-1)+s(j)),geom,chan);
    end

    %%%----- plot/save at regular intervals -----%%%
   
    if (tn > tmeasure)
        
        %%%----- store new values in array and plot -----%%%
       
        A = U(1,:);
        Au = U(2,:);
        Uz = h + B;              
       
        hovz(ihov,:) = Uz(2:end-1);
        hovh(ihov,:) = h(2:end-1);
        hovA(ihov,:) = A(2:end-1);
        hovAu(ihov,:) = Au(2:end-1);

        tmeasure = tmeasure + dtmeasure;
        ihov = ihov + 1;
              
        fp = index_fp(50);
        ct = index_city(5);
        
        fg = figure(114);
        subplot(2,2,1); 
        for k = 1:(length(s)-1)
            plot([s(k), s(k+1)],[h(k+1),h(k+1)],'b-'); hold on;
        end
        plot([s(fp), s(fp)],[0,0.04],'k:'); hold on;
        plot([s(ct), s(ct)],[0,0.04],'k:'); hold on;
        fill([chan.LR1, chan.LR2,chan.LR2,chan.LR1],[0,0,geom.hc,geom.hc],...
            'r','FaceAlpha',0.1,'LineStyle','none');
        hold off;
%         text(0.85*xmax,0.9*hmax,['t=',num2str(tn)]);
        axis([0 L 0 0.04]);
        xlabel('s','fontsize',14); ylabel('h(s,t)','fontsize',14);
        title([]);

        subplot(2,2,2);
        for k = 1:(length(s)-1)
            plot([s(k), s(k+1)],[Au(k+1),Au(k+1)],'b-'); hold on;
        end
        hold off;
        axis([0 L 0 5*max(U0(2,:))]);
        xlabel('s','fontsize',14); ylabel('Au(s,t)','fontsize',14);
        title([]);
        
        subplot(2,2,4);
        [X,Y,Xc,Yc,~] = plot_xsec_hAs(A(ct+1),s(ct),geom,chan);
        plot(Xc,Yc,'k', 'Linewidth',2); hold on;
        fill(X,Y,'b','FaceAlpha',0.1); hold off;
        text(Xc(end),0.5*geom.hr,['t=',num2str(tn)],'fontsize',14,'HorizontalAlignment', 'right');
        text(Xc(end),0.25*geom.hr,['s=',num2str(s(ct))],'fontsize',14,'HorizontalAlignment', 'right');
        
        subplot(2,2,3);
        [X,Y,Xc,Yc,hc] = plot_xsec_hAs(A(fp+1),s(fp),geom,chan);
        plot(Xc,Yc,'k', 'Linewidth',2); hold on;
        fill(X,Y,'b','FaceAlpha',0.1); hold off;
        text(Xc(end),0.5*geom.hr,['t=',num2str(tn)],'fontsize',14,'HorizontalAlignment', 'right');
        text(Xc(end),0.25*geom.hr,['s=',num2str(s(fp))],'fontsize',14,'HorizontalAlignment', 'right');
        
        set(gcf, 'Position',  [100, 100, 1000, 600])

        drawnow; pause(0.01);
        
        if (VIDEO == 1)
            frame = getframe(fg);
            writeVideo(v,frame);
        end
    end
   
end

if (VIDEO == 1)
    close(v);
end

%% SAVING DATA
if (SAVE == 1)
    
    PDE.z = hovz;
    PDE.h = hovh;
    PDE.A = hovA;
    PDE.Au = hovAu;
    
    save(fullfile(outdirs, outnames),'PDE','chan','geom');
    
end