function y = rainterm(Kkel,xx1,xx2,na,g3,nini,Lx,fac,therain,alpg,hm1,hm2)

y = 0;
for nn1 = -1:2:1 % Using simple Gaussian quadrature
    xi1    = nn1*g3;
    Na(1)  = 0.5*(1-xi1);
    Na(2)  = 0.5*(1+xi1);
    dNa(1) = -0.5;
    dNa(2) = 0.5;
    xx = xx1*Na(1)+xx2*Na(2);
    hm = hm1*Na(1)+hm2*Na(2);
    J = 0.5*Kkel;
    y = y + J*Na(na)*fac*therain-alpg*hm*dNa(na)*(hm1*dNa(1)+hm2*dNa(2))/J;
end

  
