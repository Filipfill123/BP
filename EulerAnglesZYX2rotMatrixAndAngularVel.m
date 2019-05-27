function rotMangVel = EulerAnglesZYX2rotMatrixAndAngularVel(EA)
% schema rotace ZYX (EulerAngles = [gamma;beta;alpha])
% prevod [EulerAngles,dEulerAngles,ddEulerangles] -> [rotationMatrix,omega,domega]

EulerAngles = EA(1:3,1);

gamma = EulerAngles(1);
beta = EulerAngles(2);
alpha = EulerAngles(3);

R = [cos(gamma) * cos(beta) -sin(gamma) * cos(alpha) + cos(gamma) * sin(beta) * sin(alpha) sin(gamma) * sin(alpha) + cos(gamma) * sin(beta) * cos(alpha); sin(gamma) * cos(beta) cos(gamma) * cos(alpha) + sin(gamma) * sin(beta) * sin(alpha) -cos(gamma) * sin(alpha) + sin(gamma) * sin(beta) * cos(alpha); -sin(beta) cos(beta) * sin(alpha) cos(beta) * cos(alpha);];

rotMangVel = R;

if size(EA,2) > 1
    
    dEulerAngles = EA(1:3,2);
    
    t1 = sin(gamma);
    t2 = cos(gamma);
    t3 = cos(beta);
    t4 = t2 * t3;
    t5 = t1 * t3;
    t6 = sin(beta);
    t7 = 0.1e1 / t3;
    t8 = t1 * t7;
    t7 = t2 * t7;
    H = [0 -t1 t4; 0 t2 t5; 1 0 -t6];

    omega = H*dEulerAngles;
    
    rotMangVel = [R,omega];
    
    if size(EA,2) == 3
         
        dalpha = dEulerAngles(3);
        dbeta = dEulerAngles(2);
        dgamma = dEulerAngles(1);
         
        ddEulerAngles = EA(1:3,3);
        
        dH = [0 -t2 * dgamma -t5 * dgamma - t2 * t6 * dbeta; 0 -t1 * dgamma t4 * dgamma - t1 * t6 * dbeta; 0 0 -t3 * dbeta];
        
        domega = dH*dEulerAngles + H*ddEulerAngles;
        
        rotMangVel = [R,omega,domega]; 
    end
end



