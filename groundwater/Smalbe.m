function y = Smalbe(Kkel,xx1,xx2,na,nb,g3)

y = 0;
for nn1 = -1:2:1 % Using simple Gaussian quadrature
    xi1    = nn1*g3;
    Na(1)  = 0.5*(1-xi1);
    Na(2)  = 0.5*(1+xi1);
    dNa(1) = -0.5;
    dNa(2) = 0.5;
    xx = xx1*Na(1)+xx2*Na(2);
    J = 0.5*Kkel;
    y = y+dNa(na)*dNa(nb)/J;
end
