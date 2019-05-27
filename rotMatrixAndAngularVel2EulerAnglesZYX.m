function EulerAngles = rotMatrixAndAngularVel2EulerAnglesZYX(varargin)
% schema rotace ZYX (EulerAngles = [gamma;beta;alpha])
% prevod [rotationMatrix,omega,domega] -> [EulerAngles,dEulerAngles,ddEulerangles]
% v pripade singularity (beta = +- pi/2 +- epsilon) je volen uhel alpha = 0
% (eventuelne predchozi hodnota: alpha_prev), rychlost a zrychelni volena
% pseudoinverzi matice H (Euler. kin. rce) - nebot reseni H^-1 v
% singularite neexistuje (nalezene nejlepsi mozne)

rotMangVel = varargin{1};
if nargin == 2
    alpha_prev = varargin{2};
else
    alpha_prev = 0;
end

epsilon = 1e-5;

% polohy
R = rotMangVel(1:3,1:3);

if abs(R(3,1) + 1) < epsilon % singularita
    beta = pi/2;
    alpha = alpha_prev;
    gamma = alpha_prev - atan2(R(1,2),R(1,3));
    EA = [gamma;beta;alpha];
elseif abs(R(3,1) - 1) < epsilon
    beta = -pi/2;
    alpha = alpha_prev;
    gamma = atan2(-R(1,2),-R(1,3))-alpha_prev;
    EA = [gamma;beta;alpha];
else
    sol = 1; % sol = -1 ... beta = [pi/2;3/2pi], sol = 1 ... beta = [-pi/2;pi/2] 
    beta = atan2(-R(3,1),sol*sqrt(R(3,2)^2+R(3,3)^2));
    gamma = atan2(sol*R(2,1),sol*R(1,1));
    alpha = atan2(sol*R(3,2),sol*R(3,3));
    EA = [gamma;beta;alpha];
end
EulerAngles = EA;

if size(rotMangVel,2) > 3
    % rychlosti
    omega = rotMangVel(:,4);
    
    t1 = sin(EA(1));
    t2 = cos(EA(1));
    t3 = cos(EA(2));
    t4 = t2 * t3;
    t5 = t1 * t3;
    t6 = sin(EA(2));
    t7 = 0.1e1 / t3;
    t8 = t1 * t7;
    t7 = t2 * t7;
    
    if abs(R(3,1) + 1) < epsilon % singularita
        invH = [0,0,0.5;-sin(gamma),cos(gamma),0;0,0,-0.5]; % = pinv(H)
        dEA = invH*omega;
    elseif abs(R(3,1) - 1) < epsilon
        invH = [0,0,0.5;-sin(gamma),cos(gamma),0;0,0,0.5]; % = pinv(H)
        dEA = invH*omega;
    else
        invH = [t7 * t6 t8 * t6 1; -t1 t2 0; t7 t8 0];
        dEA = invH*omega;
    end
    EulerAngles = [EA,dEA];
    
    if size(rotMangVel,2) == 5
        domega = rotMangVel(:,5);
        % zrychleni
        dH = [0 -t2 * dEA(1) -t5 * dEA(1) - t2 * t6 * dEA(2); 0 -t1 * dEA(1) t4 * dEA(1) - t1 * t6 * dEA(2); 0 0 -t3 * dEA(2)];
        ddEA = invH*(domega-dH*dEA);
        EulerAngles = [EA,dEA,ddEA];
    end   
end
   