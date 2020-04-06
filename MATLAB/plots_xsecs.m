%%% CROSS-SECTION GEOMETRY: PLOTS ETC. 

% run the first section of AuNCP_wetro0 for geom parameters

% choose s location (in city or fp)
sx = s(index_fp(5));

A1 = geom.wr*geom.hr;   % threshold area river
A2 = (geom.hr+geom.hf)*(geom.wr+geom.wf)-geom.wf*(geom.hr+0.5*geom.hf); % 2nd threshold area river
Ac = geom.wr*geom.hc; % threshold area city

h1 = geom.hr;
h2 = geom.hr + geom.hf;
hc = geom.hc;

Avec = 0:0.00001:0.004;
hvec = 0:0.0001:0.04;

dhdA = zeros(1,length(Avec));
h = zeros(1,length(Avec));
area = zeros(1,length(hvec));
W = zeros(1,length(hvec));
R = zeros(1,length(hvec));

for i = 1:length(Avec)
    [h(i), dhdA(i)] = xsec_hAs(Avec(i),sx,geom,chan);
    [area(i), W(i), R(i)] = xsec_Ahs(hvec(i),sx,geom,chan);
end


figure(100);
subplot(2,2,1);
plot(Avec, h,'k'); hold on;
if (sx > chan.LR1) && (sx < chan.LR2) % city region
    plot([Ac,Ac],[0, max(h)],'r--');
else %floodplain
    plot([A1,A1],[0, max(h)],'r--'); hold on;
    plot([A2,A2],[0, max(h)],'r--');
end
axis([Avec(1) Avec(end) 0 max(h)]);
xlabel('$A$','fontsize',14,'Interpreter','latex'); 
ylabel('$h(A,s)$','fontsize',14,'Interpreter','latex');

subplot(2,2,2);
plot(Avec, dhdA,'k'); hold on;
if (sx > chan.LR1) && (sx < chan.LR2) % city region
    plot([Ac,Ac],[0, 1.25*max(dhdA)],'r--');
else %floodplain
    plot([A1,A1],[0, 1.25*max(dhdA)],'r--'); hold on;
    plot([A2,A2],[0, 1.25*max(dhdA)],'r--');
end
axis([Avec(1) Avec(end) 0 1.25*max(dhdA)]);
xlabel('$A$','fontsize',14,'Interpreter','latex'); 
ylabel('$\partial h / \partial A$','Interpreter','latex','fontsize',14);

subplot(2,2,3);
plot(hvec, area,'k'); hold on;
if (sx > chan.LR1) && (sx < chan.LR2) % city region
    plot([hc,hc],[0, max(area)],'r--');
else %floodplain
    plot([h1,h1],[0, max(area)],'r--'); hold on;
    plot([h2,h2],[0, max(area)],'r--');
end
axis([hvec(1) hvec(end) 0 max(area)]);
xlabel('$h$','fontsize',14,'Interpreter','latex'); 
ylabel('$A(h,s)$','Interpreter','latex','fontsize',14);

subplot(2,2,4);
plot(hvec, R,'k'); hold on;
if (sx > chan.LR1) && (sx < chan.LR2) % city region
    plot([hc,hc],[0, max(R)],'r--');
else %floodplain
    plot([h1,h1],[0, max(R)],'r--'); hold on;
    plot([h2,h2],[0, max(R)],'r--');
end
axis([hvec(1) hvec(end) 0 max(R)]);
xlabel('$h$','fontsize',14,'Interpreter','latex'); 
ylabel('$R(h,s)$','fontsize',14,'Interpreter','latex');
