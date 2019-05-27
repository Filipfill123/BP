function [trans,orient,feed] = linSegPolyBlendInterp_withOrient_eval(S,data)
% na rozdil od linSegPolyBlendInterp_withOrient.m (kdy je generovan feedrate dle BAVS, AVS), zde se vraci jen zavislost na ujete draze S (oblouku) s rozlisenim res
% vypocet feedrate + vycisleni: lin seg + poly blend pro parametrizaci ujetou drahou (prirozena aproximativni/presna parametrizace)
% S ... delka ujeteho oblouku - [m]

feedType = data.feedType;
u_knot = data.u_knot; % uzlovy vektor
S_knot = data.S_knot; % ujeta draha (v translaci) v uzlovem vektoru (lin a poly segmenty)
deltaTheta = data.deltaTheta; % ujeta draha v orientaci mezi koincidencnimi body
linSegPar_trans = data.linSegPar_trans; % parametry generovane trajektorie
linSegPar_orient = data.linSegPar_orient; % parametry generovane trajektorie
polySegPar_trans = data.polySegPar_trans; % parametry generovane trajektorie
polySegPar_orient = data.polySegPar_orient; % parametry generovane trajektorie

% *** vypocet feedrate ***
S_max = S_knot(end); % max. ujeta draha po trajektorii

switch feedType
    case 1 % aproximace
        u = S/S_max;
        du_dS = 1/S_max;
        d2u_dS2 = 0;
            
    case 2 % presne
        [FEED] = linSegPolyBlend_withOrient_feedrate(S,u_knot,S_knot,polySegPar_trans);
        u = FEED(1);
        du_dS = FEED(2);
        d2u_dS2 = FEED(3);
end

% *** vycisleni trajektorie (translace + orientace) ***
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


% zavislosti na S
dp_dS = dp_du*du_dS;
d2p_dS2 = d2p_du2*du_dS^2 + dp_du*d2u_dS2;

dq_dS = dq_du*du_dS;
d2q_dS2 = d2q_du2*du_dS^2 + dq_du*d2u_dS2;

% translace: tecna rychlost
dpNorm_dS = norm(dp_dS(1:3));

% rotace: uhlova rychlost a zrychleni
rotMatVelAccel = Quaternion2rotMatrixAndAngularVel([q,dq_dS,d2q_dS2]);
R{i} = rotMatVelAccel(:,1:3);
omega = rotMatVelAccel(:,4);
domega_dS = rotMatVelAccel(:,5);

feed = [u,du_dS,d2u_dS2];
trans = [p,dp_dS,d2p_dS2];
orient = [q,dq_dS,d2q_dS2];
