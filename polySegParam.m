function [polySegPar_trans,polySegPar_orient] = polySegParam(u_knot,linSegPar_trans,linSegPar_orient,Qvect,deltaTheta)
% parametry polynomialnich segmentu
for i = 1:length(linSegPar_trans)-1
    j = 2*i;

    u1 = u_knot(j);
    u2 = u_knot(j+1);

    M = [u1^5 u1^4 u1^3 u1^2 u1 1; u2^5 u2^4 u2^3 u2^2 u2 1; 5*u1^4 4*u1^3 3*u1^2 2*u1 1 0; 5*u2^4 4*u2^3 3*u2^2 2*u2 1 0; 20*u1^3 12*u1^2 6*u1 2 0 0; 20*u2^3 12*u2^2 6*u2 2 0 0];
    I_trans = linSegPar_trans{i}(:,1)*u1 + linSegPar_trans{i}(:,2);
    O_trans = linSegPar_trans{i+1}(:,1)*u2 + linSegPar_trans{i+1}(:,2);
    dI_trans_du = linSegPar_trans{i}(:,1);
    dO_trans_du = linSegPar_trans{i+1}(:,1);
    
    theta_u1 = linSegPar_orient{i}{1}(:,1)*u1 + linSegPar_orient{i}{1}(:,2);
    theta_u2 = linSegPar_orient{i+1}{1}(:,1)*u2 + linSegPar_orient{i+1}{1}(:,2);

    uu = theta_u1/deltaTheta(i);
    duu_du = linSegPar_orient{i}{1}(:,1)/deltaTheta(i);
    d2uu_du2 = 0;
    
    SL = slerp([Qvect(:,i),Qvect(:,i+1)],uu,1);
    I_orient = SL(:,1);
    dI_orient_duu = SL(:,2);
    d2I_orient_duu2 = SL(:,3);
    dI_orient_du = dI_orient_duu*duu_du;
    d2I_orient_du2 = d2I_orient_duu2*duu_du^2;
    
    uu = theta_u2/deltaTheta(i+1);
    duu_du = linSegPar_orient{i+1}{1}(:,1)/deltaTheta(i+1);
    d2uu_du2 = 0;
    
    SL = slerp([Qvect(:,i+1),Qvect(:,i+2)],uu,1);
    O_orient = SL(:,1);
    dO_orient_duu = SL(:,2);
    d2O_orient_duu2 = SL(:,3);
    dO_orient_du = dO_orient_duu*duu_du;
    d2O_orient_du2 = d2O_orient_duu2*duu_du^2;

    for k = 1:size(O_trans,1)
        b_trans = [I_trans(k);O_trans(k);dI_trans_du(k);dO_trans_du(k);0;0];
        polySegPar_trans{i}(k,:) = (M^-1*b_trans)';
    end
    
    for k = 1:size(O_orient,1)
        b_orient = [I_orient(k);O_orient(k);dI_orient_du(k);dO_orient_du(k);d2I_orient_du2(k);d2O_orient_du2(k)];
        polySegPar_orient{i}(k,:) = (M^-1*b_orient)';
    end
end