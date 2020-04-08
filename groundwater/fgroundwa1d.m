%
%
% Groundwater FEM code
%
% M_ij h_j^n+1+(Lc/m_por*sige)*h_cm^n+1 =
% M_ij h_j^n+(Lc/m_por*sige)*h_cm^n + dt*b_i^n-(dt/m_por*sige)*sqrt(g)*max(2*h^n_cm,0)^3/2
% b_i^n = int_0^Ly (alpha*g*h^n_m d_y varphi d_y h^n_m) dy
%
% FEM code with Dirichlet boundary conditions
%
% parameters:
ntimstepper = 0;		   
g = 9.81; % m/s
Ly = 0.85; % m
Lc = -0.05; % m
Lca = abs(Lc);
wv = 0.1; % width Hele-Shaw moor cell in m
sigm = 0.8;    % fraction of moor pores filled
sigm0 = 0.1;   % fraction of water remaining in moor pores after water has seeped out
sigme = sigm;  % 
mpor = 0.3;    % porosity moor
nu = 10^(-6);  % viscosity water m^2/s
kperm = 10^-8; % permeability
alph = kperm/(mpor*nu*sigm);
Rain0 = 0.000125; % m/s
Rain = Rain0;
%
Nk = 20;   % number of elements
Nn = Nk+1; % number of nodes
Mm = zeros(Nn,Nn);
Sm = zeros(Nn,Nn);
ngrid = 1;
switch ngrid
case 1 % regular grid
  dy = Ly/Nk;
  yy = transpose(0:dy:Ly);  % Regular FEM grid
end
Kk = (yy(2:Nn)-yy(1:Nn-1)); % finite element lengths
%
%
%
TimeUnit = 10; % 10; %s
Vrate = Ly*wv*Rain0*[1,2,4,8,18]*TimeUnit*1000; % flow rates in liters/TimeUnit
Tend = 15*TimeUnit; % 15*TimeUnit;
Nt = 5;
dtmeas = TimeUnit/Nt;
tmeas = dtmeas;
ntimestepper = 1; 
g3 = sqrt(1.0/3.0);
nini = 3;
%
% Make mass and Laplace matrices using matrix assembly as wella s initial condition
%
switch nini
case 1 %
 hcm = 0.0;
case 2 %
 hcm = 0.0;
case 3 %
  hcm = 0.0;
end
%
bini = zeros(Nn,1);
for nk=1:Nk
   for na = 1:2
       ii = (nk-1)+na;  
       bini(ii) = bini(ii)+bgina(Kk(nk),yy(nk),yy(nk+1),na,g3,nini,Ly); 
       for nb = 1:2
           jj = (nk-1)+nb;
           Mm(ii,jj) = Mm(ii,jj)+Mmalbe(Kk(nk),yy(nk),yy(nk+1),na,nb,g3);
           Sm(ii,jj) = Sm(ii,jj)+Smalbe(Kk(nk),yy(nk),yy(nk+1),na,nb,g3);
      end	 		    
   end
end
%
% Initial condition and plotting
%
timee = 0.0;
uu = zeros(Nn,1);
uu(1:Nn) = Mm\bini;
%
%
figno = 21;
figure(figno); clf;
subplot(2,1,1);
plot(yy,uu,'k','linewidth',2);
xlabel('y','fontsize',18);
ylabel('h_m(y,t)','fontsize',18);
subplot(2,1,2);
plot(timee,hcm,'o' ); hold on;
xlabel('t','fontsize',18);
ylabel('h_{cm}(t)','fontsize',18);
drawnow;
% umin = min(uu);
% umax = max(uu);
% axis([0 Ly umin umax]);
pause(2);
%
% Time stepping
%
nrain = 1;
%
ntimstepper = 1;
switch ntimestepper
case 1
  CFL = 0.3;
  dt = CFL*0.5*dy^2; % FD estimate 
case 2
  CFL = 10;
  dt = CFL*0.5*dy^2; % FD estimate 
end
%
% tmeas = dt;
% dtmeas = dt;
%
Ml = zeros(Nn,Nn);
Ml = Mm;
Ml(1,1) = Ml(1,1)+Lca/(mpor*sigme);
while (timee <= Tend)
    %
    timee = timee + dt;
    % 
    bini = zeros(Nn,1); % size(bini)
    for nk=1:Nk
       for na = 1:2
          ii = (nk-1)+na;  
          bini(ii) = bini(ii)+rainterm(Kk(nk),yy(nk),yy(nk+1),na,g3,nini,Ly,1.0/(mpor*sigme),nrain*Rain0,alph*g,uu(nk),uu(nk+1)); 
        end
    end
    hcmold = 1.0;
    % size(bini)
    % switch ntimestepper
    % case 1 % forward Euler
      %% unew(2:Nn) = Mm(2:Nn,2:Nn)\( Mm(2:Nn,2:Nn)*uu(2:Nn)+dt*bini(2:Nn) );
      % faulty hcmn = ( hcm - (dt/Lca)*sqrt(g)*max(2*hcm/3,0.0)^(1.5) )/(1.0 - dt*(0.5*alph*mpor*sigme*g)/(Lca)*0.5*Kk(1)*(uu(2)-hcm) );
      %% hcmn = hcm - (dt/Lca)*sqrt(g)*max(2*hcm/3,0.0)^(1.5) + dt*(0.5*alph*mpor*sigme*g*hcm)/(Lca)*0.5*Kk(1)*(uu(2)-hcm) 
      %bini(1) = -(dt/Lca)*sqrt(g)*max(2*hcm/3,0.0)^(1.5);
      %Ml(1,1) = 1+dt*(0.5*alph*mpor*sigme*g*hcm)/(Lca)*0.5*Kk(1);
      %Ml(1,2) = -dt*(0.5*alph*mpor*sigme*g*hcm)/(Lca)*0.5*Kk(1);
      bini(1) = bini(1)-sqrt(g)*max(2*hcm/3,0.0)^(1.5)/(mpor*sigme);
      unew = Ml\( Ml*uu+dt*bini );
    % end
    % pause;
    hcm = unew(1);
    uu(1) = hcm;
    uu(2:Nn) = transpose(unew(2:Nn));
    % plotting
    if ( timee > tmeas )
      timee
      tmeas = tmeas + dtmeas;
      figure(figno); % clf;
      subplot(2,1,1);
      plot(yy,uu,'k','linewidth',2);
      xlabel('y','fontsize',18); 
      ylabel('h_m(y,t)','fontsize',18);
      % axis([0 Ly umin umax]);
      subplot(2,1,2);
      plot(timee,hcm,'o' ); hold on;
      drawnow;
      pause(0.2);
    end
    %
end
