function [trans,orient] = linSegPolyBlendEval_withOrient(u,u_knot,linSegPar_trans,linSegPar_orient,polySegPar_trans,polySegPar_orient,deltaTheta)
% vycisleni: lin seg + poly blend

seg = findSpan(u,u_knot);

if mod(seg,2) == 0
    % poly blend
    i = seg/2;
    
    % translace
    p = polySegPar_trans{i}(:,1)*u^5 + polySegPar_trans{i}(:,2)*u^4 + polySegPar_trans{i}(:,3)*u^3 + polySegPar_trans{i}(:,4)*u^2 + polySegPar_trans{i}(:,5)*u + polySegPar_trans{i}(:,6);
    dp_du = 5*polySegPar_trans{i}(:,1)*u^4 + 4*polySegPar_trans{i}(:,2)*u^3 + 3*polySegPar_trans{i}(:,3)*u^2 + 2*polySegPar_trans{i}(:,4)*u + polySegPar_trans{i}(:,5);
    d2p_du2 = 20*polySegPar_trans{i}(:,1)*u^3 + 12*polySegPar_trans{i}(:,2)*u^2 + 6*polySegPar_trans{i}(:,3)*u + 2*polySegPar_trans{i}(:,4);
    
    % orientace (kvaternion)
    q_unnorm = polySegPar_orient{i}(:,1)*u^5 + polySegPar_orient{i}(:,2)*u^4 + polySegPar_orient{i}(:,3)*u^3 + polySegPar_orient{i}(:,4)*u^2 + polySegPar_orient{i}(:,5)*u + polySegPar_orient{i}(:,6);
    dq_du_unnorm = 5*polySegPar_orient{i}(:,1)*u^4 + 4*polySegPar_orient{i}(:,2)*u^3 + 3*polySegPar_orient{i}(:,3)*u^2 + 2*polySegPar_orient{i}(:,4)*u + polySegPar_orient{i}(:,5);
    d2q_du2_unnorm = 20*polySegPar_orient{i}(:,1)*u^3 + 12*polySegPar_orient{i}(:,2)*u^2 + 6*polySegPar_orient{i}(:,3)*u + 2*polySegPar_orient{i}(:,4);    
    
    % normovani
    norm_q_unnorm = norm(q_unnorm);
    dnorm_q_unnorm_du = 1/norm_q_unnorm*q_unnorm'*dq_du_unnorm;
    d2norm_q_unnorm_du2 = -1/norm_q_unnorm^2*dnorm_q_unnorm_du*q_unnorm'*dq_du_unnorm + 1/norm_q_unnorm*dq_du_unnorm'*dq_du_unnorm + 1/norm_q_unnorm*q_unnorm'*d2q_du2_unnorm;

    q = q_unnorm/norm_q_unnorm;
    dq_du = (dq_du_unnorm*norm_q_unnorm - q_unnorm*dnorm_q_unnorm_du)/(norm_q_unnorm^2);
    d2q_du2 = ((d2q_du2_unnorm*norm_q_unnorm + dq_du_unnorm*dnorm_q_unnorm_du - dq_du_unnorm*dnorm_q_unnorm_du - q_unnorm*d2norm_q_unnorm_du2)*norm_q_unnorm^2 - 2*(dq_du_unnorm*norm_q_unnorm - q_unnorm*dnorm_q_unnorm_du)*q_unnorm'*dq_du_unnorm)/(norm_q_unnorm^4);
else
    % lin. segment
    i = floor(seg/2)+1;
    
    % translace
    p = linSegPar_trans{i}(:,1)*u + linSegPar_trans{i}(:,2);
    dp_du = linSegPar_trans{i}(:,1);
    d2p_du2 = 0*dp_du;
    
    % orientace
    theta = linSegPar_orient{i}{1}(:,1)*u + linSegPar_orient{i}{1}(:,2);
    uu = theta/deltaTheta(i);
    duu_du = linSegPar_orient{i}{1}(:,1)/deltaTheta(i);
    d2uu_du2 = 0;
    
    SL = slerp([linSegPar_orient{i}{2}(:,1),linSegPar_orient{i}{2}(:,2)],uu,1); % slerp([Quat1,Quat2],uu), uu = <0,1>
    q = SL(:,1);
    dq_duu = SL(:,2);
    d2q_duu2 = SL(:,3);  
    dq_du = dq_duu*duu_du;
    d2q_du2 = d2q_duu2*duu_du^2;
end 

trans = [p,dp_du,d2p_du2];
orient = [q,dq_du,d2q_du2];