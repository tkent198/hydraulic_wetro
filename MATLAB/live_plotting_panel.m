% NOTE: legends are SLOW for live plotting

fs = 14; %fontsize

% Panel: 1 = 2x5; 2 = 3x4
panel = 2;

switch panel 
    case 1
        % cols rows
        nrow = 2;
        ncol = 5;
        %placers
        pr = 1;
        phist = 2;
        pm = 3; 
        pc = 4;
        pQh = 6;
        pht = 7;
        phct = 9;
        phfp = 8;
        pAu = 5;
        ph = 10;
        %subaxis
        sp = 0.05;
        pa = 0.;
        ma = 0.075;
    case 2
        % cols rows
        nrow = 3;
        ncol = 4;
        %placers
        pr = 1;
        phist = 2;
        pm = 5; 
        pc = 6;
        pQh = 9;
        pht = 10;
        phct = 12;
        phfp = 11;
        pAu = [3 4];
        ph = [7 8];
        %subaxis
        sp = 0.06;
        pa = 0.;
        ma = 0.05;
end

% mainfig = figure(111);
% mainfig.WindowState = 'maximized';

% subplot(2,2,1);
% subplot(2,6,2);
% subplot(2,5,3);
subaxis(nrow,ncol,pc, 'Spacing', sp, 'Padding', pa, 'Margin', ma);
plot(tijd,h1c,'ok','linewidth',1); hold on;
plot(tijd,h2c,'*c','linewidth',1); hold on;
plot(tijd,h3c,'xb','linewidth',1); hold on;
plot(tijd,hres/10,'*r','linewidth',1); hold on;
legend({'Canal 1','Canal 2', 'Canal 3', 'Res/10'},'Location','northwest','fontsize',fs);
ylim([0,0.025]);
xlabel('t [s]','fontsize',fs);
ylabel('Water depth h [m]', 'fontsize',fs); hold on;

% subplot(2,2,2);
% subplot(2,6,4);
% subplot(2,5,8);
subaxis(nrow,ncol,pht, 'Spacing', sp, 'Padding', pa, 'Margin', ma);
plot([tijdo,tijd],[hro,h(1)],'b','linewidth',1); hold on; %s=0
plot([tijdo,tijd],[hrLo,h(nxLc1+2)],'k','linewidth',1); hold on; %s=L_c1
plot([tijdo,tijd],[0.02,0.02],':r','linewidth',1); hold on; %h flood?
% axis([0 Tend 0 max(h2c)]);
% legend({'h(s=0,t)','h(s=L_{1c},t)'},'Location','southeast','fontsize',fs);
ylim([0,0.025]);
xlabel('t [s]', 'fontsize',fs);
ylabel('River level h [m]', 'fontsize',fs); hold on;

% subplot(2,2,3);
% subplot(2,6,3);
% subplot(2,5,2);
subaxis(nrow,ncol,pm, 'Spacing', sp, 'Padding', pa, 'Margin', ma);
plot(yy,hm,'-','linewidth',2);
% legend({'Groundwater level'},'Location','southeast','fontsize',fs);
ylim([0,0.12]);
xlabel('y [m]', 'fontsize',fs);
ylabel('h_m(y,t) [m]', 'fontsize',fs);

% subplot(2,2,4);
% subplot(2,6,1);
% subplot(2,5,1);
subaxis(nrow,ncol,pr, 'Spacing', sp, 'Padding', pa, 'Margin', ma);
plot(tijd,Rm(1)/Rain0,'oc','linewidth',1); hold on;
plot(tijd,Rr(1)/Rain0,'db','linewidth',1); hold on;
plot(tijd,(Rm(1)+Rr(1))/Rain0,'xr','linewidth',1); hold  on;
legend({'Moor','Res', 'Both'},'Location','northeast','fontsize',fs);
ylim([0,18]);
xlabel('t [s]', 'fontsize',fs);
ylabel('Rainfall factor', 'fontsize',fs);
%

% fg = figure(114);
% subplot(2,2,1);
% subplot(2,6,[11 12]);
% subplot(2,5,9);
subaxis(nrow,ncol,ph, 'Spacing', sp, 'Padding', pa, 'Margin', ma);
for k = 1:(length(s)-1)
    plot([s(k), s(k+1)],[h(k+1),h(k+1)],'b-', 'Linewidth',1); hold on;
end
plot([chan.s_r, chan.s_r],[0,0.04],'r:'); hold on;
plot([chan.s_m, chan.s_m],[0,0.04],'r:'); hold on;
plot([Lc1, Lc1],[0,0.04],'r:'); hold on;
plot([s(fp), s(fp)],[0,0.04],'k:'); hold on;
plot([s(ct), s(ct)],[0,0.04],'k:'); hold on;
fill([chan.LR1, chan.LR2,chan.LR2,chan.LR1],[0,0,geom.hc,geom.hc],...
    'r','FaceAlpha',0.1,'LineStyle','none');
hold off;
% text(0.85*xmax,0.9*hmax,['t=',num2str(tn)]);
axis([0 L 0.01 0.03]);
xlabel('s','fontsize',fs); ylabel('h(s,t)','fontsize',fs);
title([]);

% subplot(2,2,2);
% subplot(2,6,[5 6]);
% subplot(2,5,10);
subaxis(nrow,ncol,pAu, 'Spacing', sp, 'Padding', pa, 'Margin', ma);
for k = 1:(length(s)-1)
    plot([s(k), s(k+1)],[Au(k+1),Au(k+1)],'b-', 'Linewidth',1); hold on;
end
plot([chan.s_r, chan.s_r],[0,0.04],'r:'); hold on;
plot([chan.s_m, chan.s_m],[0,0.04],'r:'); hold on;
plot([Lc1, Lc1],[0,0.04],'r:'); hold on;
hold off;
axis([0 L 0.0001 0.0004]);
xlabel('s','fontsize',fs); ylabel('Au(s,t)','fontsize',fs);
title([]);

% subplot(2,2,4);
% subplot(2,6,10);
% subplot(2,5,7);
subaxis(nrow,ncol,phct, 'Spacing', sp, 'Padding', pa, 'Margin', ma);
[X,Y,Xc,Yc,~] = plot_xsec_hAs(A(ct+1),s(ct),geom,chan);
plot(Xc,Yc,'k', 'Linewidth',2); hold on;
fill(X,Y,'b','FaceAlpha',0.1); hold off;
text(Xc(1)+0.01,0.9*Yc(1),['t=',num2str(tijd)],'fontsize',fs,'HorizontalAlignment', 'left');
text(Xc(1)+0.01,0.8*Yc(1),['s=',num2str(s(ct))],'fontsize',fs,'HorizontalAlignment', 'left');
axis([min(Xc) max(Xc) min(Yc) max(Yc)]);

% subplot(2,2,3);
% subplot(2,6,9);
% subplot(2,5,6);
subaxis(nrow,ncol,phfp, 'Spacing', sp, 'Padding', pa, 'Margin', ma);
[X,Y,Xc,Yc,hc] = plot_xsec_hAs(A(fp+1),s(fp),geom,chan);
plot(Xc,Yc,'k', 'Linewidth',2); hold on;
fill(X,Y,'b','FaceAlpha',0.1); hold off;
text(Xc(1)+0.01,0.9*Yc(1),['t=',num2str(tijd)],'fontsize',fs,'HorizontalAlignment', 'left');
text(Xc(1)+0.01,0.8*Yc(1),['s=',num2str(s(fp))],'fontsize',fs,'HorizontalAlignment', 'left');
axis([min(Xc) max(Xc) min(Yc) max(Yc)]);

%         figure(13);
% subplot(1,2,1);
% subplot(2,6,7);
% subplot(2,5,4);
subaxis(nrow,ncol,phist, 'Spacing', sp, 'Padding', pa, 'Margin', ma);
% ntot = Tend/TimeUnit;
histogram(nchisto, 'Normalization','pdf'); hold on;
plot(rainfac,rainpdf,'hk','Linewidth',2); hold off;
xlabel('Rain / wd', 'fontsize',fs);
ylabel('Density', 'fontsize',fs);
axis([-0.5 20 0 0.4]);

% subplot(2,6,8);
% subplot(2,5,5);
subaxis(nrow,ncol,pQh, 'Spacing', sp, 'Padding', pa, 'Margin', ma);
plot(h(nxLc1+2), U(2,nxLc1+2),'xk'); hold on;
text(0.011,0.000275,['s=',num2str(s(ct))],'fontsize',fs,'HorizontalAlignment', 'left'); 
xlabel('h', 'fontsize',fs);
ylabel('Q', 'fontsize',fs);
axis([0.01 0.03 0.0001 0.0003]);

drawnow; pause(0.2);

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