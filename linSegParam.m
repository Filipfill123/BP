function [linSegPar_trans,linSegPar_orient] = linSegParam(u_knot,Pvect,Ovect,Ivect,deltaThetaOvect,deltaThetaIvect,deltaTheta,Qvect)
% parametry linearnich segmentu
% TRANS pro i-ty segment: [x,y,z] = Lp1_trans*u + Lp2_trans, Lp_trans{i} = [Lp1_trans,Lp2_trans]
% ORIENT pro i-ty segment theta = Lp1_orient*u + Lp2_orient, Lp_orient{i} = [Lp1_orient,Lp2_orient], theta lezi v int <0,deltaTheta(i)> a odpovida pocatku a konci lin. useku pri pouziti SLERP interpolace
numOfPoints = size(Pvect,2);

for i = 1:numOfPoints - 1
    if i == 1
        j = 2*(i-1)+1;
        Lp1_trans = (Ivect(:,i) - Pvect(:,i))/(u_knot(j+1)-u_knot(j));
        Lp2_trans = Pvect(:,i) - (Ivect(:,i) - Pvect(:,i))/(u_knot(j+1)-u_knot(j))*u_knot(j);
        
        Lp1_orient = (deltaThetaIvect(:,i) - 0)/(u_knot(j+1)-u_knot(j));
        Lp2_orient = 0 - (deltaThetaIvect(:,i) - 0)/(u_knot(j+1)-u_knot(j))*u_knot(j);
    elseif i == numOfPoints-1
        j = 2*(i-1)+1;
        Lp1_trans = (Pvect(:,i+1) - Ovect(:,i-1))/(u_knot(j+1)-u_knot(j));
        Lp2_trans = Ovect(:,i-1) - (Pvect(:,i+1) - Ovect(:,i-1))/(u_knot(j+1)-u_knot(j))*u_knot(j);
        
        Lp1_orient = (deltaTheta(:,i) - deltaThetaOvect(:,i-1))/(u_knot(j+1)-u_knot(j));
        Lp2_orient = deltaThetaOvect(:,i-1) - (deltaTheta(:,i) - deltaThetaOvect(:,i-1))/(u_knot(j+1)-u_knot(j))*u_knot(j);
    else
        j = 2*(i-1)+1;
        Lp1_trans = (Ivect(:,i) - Ovect(:,i-1))/(u_knot(j+1)-u_knot(j));
        Lp2_trans = Ovect(:,i-1) - (Ivect(:,i) - Ovect(:,i-1))/(u_knot(j+1)-u_knot(j))*u_knot(j);
        
        Lp1_orient = (deltaThetaIvect(:,i) - deltaThetaOvect(:,i-1))/(u_knot(j+1)-u_knot(j));
        Lp2_orient = deltaThetaOvect(:,i-1) - (deltaThetaIvect(:,i) - deltaThetaOvect(:,i-1))/(u_knot(j+1)-u_knot(j))*u_knot(j);
    end

    linSegPar_trans{i} = [Lp1_trans,Lp2_trans];
    linSegPar_orient{i} = {[Lp1_orient,Lp2_orient],[Qvect(:,i),Qvect(:,i+1)]};
end
