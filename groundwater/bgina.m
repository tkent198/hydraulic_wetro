function y = bgina(Kkel,xx1,xx2,na,g3,nini,Lx)

y = 0;
for nn1 = -1:2:1 % Using simple Gaussian quadrature
    xi1    = nn1*g3;
    Na(1)  = 0.5*(1-xi1);
    Na(2)  = 0.5*(1+xi1);
    xx = xx1*Na(1)+xx2*Na(2);
    switch nini
    case 1 % 
       u0 = 0.0;
    case 2 
       u0 = 0.0;
    case 3
       u0 = 0.0;
    end
    J = 0.5*Kkel;
    y = y+J*Na(na)*u0;
end

  
