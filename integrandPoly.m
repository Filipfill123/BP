function [integrand,dintegrand_du] = integrandPoly(u,polyPar)
% integrand pro dany segment polyPar (normalized - u = <0,1>!!)
% pouze z prvnich 3 souradnic (x,y,z)
for j = 1:length(u)
    dp_du = 5*polyPar(:,1)*u(j)^4 + 4*polyPar(:,2)*u(j)^3 + 3*polyPar(:,3)*u(j)^2 + 2*polyPar(:,4)*u(j) + polyPar(:,5);
    d2p_du2 = 20*polyPar(:,1)*u(j)^3 + 12*polyPar(:,2)*u(j)^2 + 6*polyPar(:,3)*u(j) + 2*polyPar(:,4);
    integrand(j) = norm(dp_du(1:3));
    dintegrand_du(j) = dp_du(1:3)'/integrand(j)*d2p_du2(1:3);
end     
