function varargout = forwardKinematics(varargin)
% vypocet prime kinematiky serioveho manipulatoru 
% 
% 
% varargin = {jointCoords,kinPar,N}
%   jointCoords = [q,dq,ddq]
% 
%   kinPar ... kinematicke parametry manipulatoru
%   kinPar.DHpar = [d_1,theta_1,a_1,alpha_1; d_2,theta_2,a_2,alpha_2;... ]
%   kinPar.jointType = [joint1, joint2,...] = ['P','R',...]
%   kinPar.Qhome = [q1h; q2h; ...]
%   kinPar.baseComp = [O0b,R0b] ... posunuti a orientace (matice rotace) -
%       kompenzace polohy zakladny
%   kinPar.endEffComp = [OeN,ReN] ... posunuti a orientace (matice rotace) -
%       kompenzace polohy koncoveho efektoru
% 
%   pro N = i  ... poloha  i-teho s.s. (kompenzace polohy konc. eff. "endEffComp" je ignorovana!)
%   pokud N neni zadano => poloha s.s. konc. efektoru Fe
% 
% varargout = {genCoords}
%   genCoords = [X,dX,ddX]
%   obecny a defaultni vystup je uvazovan jako:  
%   X = [O_n^0,R_n^0]
%   dX = [dO_n^0,omega_n^0]
%   ddX = [ddO_n^0,domega_n^0]
%   - muze byt UZIVATELEM zmeneno na neco jineho - jinak definovane
%       zobecnene souradnice

jointCoords = varargin{1};
kinPar = varargin{2};
DHpar = kinPar.DHpar;
jointType = kinPar.jointType;
Qhome = kinPar.Qhome;
endEffComp = kinPar.endEffComp;
baseComp = kinPar.baseComp;

% kompenzace polohy konc. ef.
Oen = endEffComp(:,1);
Ren = endEffComp(:,2:4);
TeN = [[Ren,Oen];[0,0,0,1]];

% kompenzace polohy zakladny
O0b = baseComp(:,1);
R0b = baseComp(:,2:4);
T0b = [[R0b,O0b];[0,0,0,1]];

if nargin > 2 % delka kin. retezce
    N = abs(varargin{3});
    endEffComp_active = false;
else
    N = size(DHpar,1);
    endEffComp_active = true;
end

% POLOHY
% vypocet hom. transformacnich matic
q = jointCoords(1:N,1);

% kompenzace polohy zakladny
TT{1} = T0b; %TT{i} = T_{i-1}^0 => TT = {T0b,T1b,T2b,T3b,...,TNb,Teb}

for i = 1:N
    if strcmp(jointType(i),'P')
        DHpar(i,1) = q(i);
    else
        DHpar(i,2) = q(i);
    end
    
    TT{i+1} = TT{i}*DH(DHpar(i,:));
end

% kompenzace polohy konc. efektoru
if endEffComp_active
    TT{N+2} = TT{N+1}*TeN;
end

X = [TT{end}(1:3,4),TT{end}(1:3,1:3)]; % [O,R]

% RYCHLOSTI, ZRYCHLENI
dq = jointCoords(1:N,2);
ddq = jointCoords(1:N,3);

if endEffComp_active
    [J,dJ,Jcomp,Jadd] = kinJacobian([q,dq],DHpar(1:N,:),jointType(1:N),baseComp,endEffComp);
else
    [J,dJ,Jcomp,Jadd] = kinJacobian([q,dq],DHpar(1:N,:),jointType(1:N),baseComp,[zeros(3,1),eye(3,3)]);
end

dX = Jcomp*J*dq; % [dO;omega]
dX = [dX(1:3),dX(4:6)]; % [dO,omega]

ddX = Jcomp*(dJ*dq + J*ddq) + Jadd; % [ddO;domega]
ddX = [ddX(1:3),ddX(4:6)]; % [ddO,domega]

% prepocet na X_,dX_,ddX_, X_ = [x;y;phi]
X_ = [X(1:2,1);sum(q)];
dX_ = [dX(1:2,1);sum(dq)];
ddX_ = [ddX(1:2,1);sum(ddq)];

varargout{1} = [X_,dX_,ddX_];

