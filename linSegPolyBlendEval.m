function [p,dp_du,d2p_du2] = linSegPolyBlendEval(u,u_knot,linSegPar,polySegPar)
% vycisleni: lin seg + poly blend

seg = findSpan(u,u_knot);

if mod(seg,2) == 0
    % poly blend
    i = seg/2;
    p = polySegPar{i}(:,1)*u^5 + polySegPar{i}(:,2)*u^4 + polySegPar{i}(:,3)*u^3 + polySegPar{i}(:,4)*u^2 + polySegPar{i}(:,5)*u + polySegPar{i}(:,6);
    dp_du = 5*polySegPar{i}(:,1)*u^4 + 4*polySegPar{i}(:,2)*u^3 + 3*polySegPar{i}(:,3)*u^2 + 2*polySegPar{i}(:,4)*u + polySegPar{i}(:,5);
    d2p_du2 = 20*polySegPar{i}(:,1)*u^3 + 12*polySegPar{i}(:,2)*u^2 + 6*polySegPar{i}(:,3)*u + 2*polySegPar{i}(:,4);
else
    % lin. segment
    i = floor(seg/2)+1;
    p = linSegPar{i}(:,1)*u + linSegPar{i}(:,2);
    dp_du = linSegPar{i}(:,1);
    d2p_du2 = 0*dp_du;
end 
