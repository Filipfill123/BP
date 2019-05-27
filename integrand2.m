function in = integrand2(u,u_knot,linSegPar,polySegPar)

for i = 1:length(u)
    [p_int,dp_int_du,d2p_int_du2] = linSegPolyBlendEval(u(i),u_knot,linSegPar,polySegPar);
    in(i) = norm(dp_int_du);
end