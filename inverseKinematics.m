function varargout = inverseKinematics(varargin)
% vypocet inverzni kinematiky serioveho manipulatoru 
% IMPLEMETUJE UZIVATEL - k vypoctu zavislosti rychlosti,zrychlenich vhodne vyuzit 
%   kinJacobian.m
% 
% varargin = {genCoords,kinPar}
%   kinPar ... kinematicke parametry manipulatoru
%   kinPar.DHpar = [d_1,theta_1,a_1,alpha_1; d_2,theta_2,a_2,alpha_2;... ]
%   kinPar.jointType = [joint1, joint2,...] = ['P','R',...]
%   kinPar.Qhome = [q1h; q2h; ...]
%   kinPar.baseComp = [O0b,R0b] ... posunuti a orientace (matice rotace) -
%       kompenzace polohy zakladny
%   kinPar.endEffComp = [OeN,ReN] ... posunuti a orientace (matice rotace) -
%       kompenzace polohy koncoveho efektoru
% 
%   genCoords = [X,dX,ddX]
%   X ... zobecnene souradnice (specifikuje uzivatel - napr. [O_N^0,R_N^0],...)
% 
% varargout = {q,dq,ddq}
% 
genCoords = varargin{1}; % [X,dX,ddX], X = [x;y;phi]
kinPar = varargin{2};
sol = varargin{3}; % sol = -1 (S2 < 0), sol = 1 (S2 > 0)

DHpar = kinPar.DHpar;
jointType = kinPar.jointType;
Qhome = kinPar.Qhome;

L1 = DHpar(1,3);
L2 = DHpar(2,3);
L3 = DHpar(3,3);

% polohy
x = genCoords(1,1);
y = genCoords(2,1);
phi = genCoords(3,1);

Sphi = sin(phi);
Cphi = cos(phi);
C2 = ((x - L3*Cphi)^2 + (y - L3*Sphi)^2 - L1^2 - L2^2)/(2*L1*L2);
S2 = sol*realsqrt(1-C2^2);
theta2 = atan2(S2,C2);

S1 = -(L1 * L3 * Sphi - L1 * y + C2 * L2 * L3 * Sphi - C2 * y * L2 - L3 * Cphi * S2 * L2 + L2 * x * S2) / (L2 ^ 2 + 0.2e1 * L2 * L1 * C2 + L1 ^ 2);
C1 = -(L3 * Sphi * S2 * L2 - y * S2 * L2 + L1 * L3 * Cphi - x * L1 + C2 * L2 * L3 * Cphi - x * C2 * L2) / (L2 ^ 2 + 0.2e1 * L2 * L1 * C2 + L1 ^ 2);
theta1 = atan2(S1,C1);

theta3 = phi - theta1 - theta2;

Q = [theta1;theta2;theta3];

% rychlosti
dX = genCoords(:,2); %dX = [dx;dy;dphi]
[J] = kinJacobian(Q,DHpar,jointType); % [dO,omega] = J*dQ, dO = [dx;dy], omega = [0;0;sum(dQ)]
Jred = [J(1:2,:);[1,1,1]]; % dX = Jred*dQ
dQ = Jred^-1*dX;

% zrychleni
ddX = genCoords(:,3); %ddX = [ddx;ddy;ddphi]
[J,dJ] = kinJacobian([Q,dQ],DHpar,jointType); % [ddO,domega] = dJ*dQ + J*ddQ , ddO = [ddx;ddy], domega = [0;0;sum(ddQ)]
Jred = [J(1:2,:);[1,1,1]]; 
dJred = [dJ(1:2,:);[0,0,0]]; % ddX = dJred*dQ + Jred*ddQ
ddQ = Jred^-1*(ddX - dJred*dQ);

varargout{1} = [Q,dQ,ddQ];