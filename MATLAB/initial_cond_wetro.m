function  [U0, B, h0] = initial_cond_wetro(s,Neq,Nk,geom,chan)

%%% Initial conditions for wetro0 system

% INPUT:
% s = cell edge coordinates
% Neq = no. of equations in the system
% Nk = no. of cells
% geom  = geometry pars
% chan = domain pars

% OUTPUT:
% U0: [Neq, Nk] array consisting of FV [A,Au] data
% B: [1, Nk] array consisting of FV topography data
% h0: [1, Nk] array consisting of FV h data


%%%---------- ICs for the Neq equations ----------------%%%

%%% for depth h: flat
ic_h = 0.0135*ones(1,length(s));
% ic_h = 0.5*(2*geom.hr+geom.hf)*ones(1,length(s));
% ic_h = 1.2*(geom.hr+geom.hf)*ones(1,length(s));

%%% for area A:
ic_A = zeros(1,length(s));
W = zeros(1,length(s));
R = zeros(1,length(s));

for i = 1:Nk+1
    [ic_A(i), W(i), R(i)] = xsec_Ahs(ic_h(i),s(i),geom,chan);
end

%%% kinetic solution for u(s,0)
ic_u = R.^(2/3)*sqrt(-geom.dbds)/geom.Cm;
% ic_u = ic_u*ones(1,length(s));


%%% for b: linear with gradient dbds
L = chan.LR3;
bc = -L*geom.dbds;
B = linspace(bc,0,Nk+1);


% Define array and fill with FV (piecewise constant) initial data 
U0 = zeros(Neq,Nk);
B = 0.5*(B(1:Nk) + B(2:Nk+1)); % b
h0 = 0.5*(ic_h(1:Nk) + ic_h(2:Nk+1)); % h
U0(1,:) = 0.5*(ic_A(1:Nk)+ic_A(2:Nk+1)); % A
U0(2,:) = 0.5*(ic_A(1:Nk).*ic_u(1:Nk) + ic_A(2:Nk+1).*ic_u(2:Nk+1)); % Au


end