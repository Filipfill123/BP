function varargout = inverseDynamicModel(varargin)
% vypocita inverzni dynamicky model manipulatoru
%
% varargin = [jointCoords,endEff_force_moment,kinPar,dynPar]
%
% jointCoords = [q,dq,ddq] ... polohy, rychlosti, zrychleni kl. souradnic
% endEff_force_moment = [f;u] ... sila, moment pusobici na konc. eff. (s.s.
%   Fe) w.r.t. s.s. Fb
% kinPar
%   kinPar.DHpar = [d_1,theta_1,a_1,alpha_1; d_2,theta_2,a_2,alpha_2;... ]
%   kinPar.jointType = [joint1, joint2,...] = ['P','R',...]
%   kinPar.Qhome = [q1h; q2h; ...]
%   kinPar.baseComp = [O0b,R0b] ... posunuti a orientace (matice rotace) -
%       kompenzace polohy zakladny
%   kinPar.endEffComp = [OeN,ReN] ... posunuti a orientace (matice rotace) -
%       kompenzace polohy koncoveho efektoru
% dynPar
%   dynPar.mass = [m1; m2; ...]
%   dynPar.inertiaTensor = [I1, I2, ...] ... w.r.t. s.s. prislusneho 
%       ramene (F1,F2,...) => konstantni (nezavisi na poloze manip.!)
%   dynPar.gravityCenter = [T1, T2, ...] ... w.r.t. s.s. prislusneho 
%       ramene (F1,F2,...) => konstantni (nezavisi na poloze manip.!)
%   dynPar.gravityVector = [0;0;-9.81];
%
% varargout = [tau1;tau2;...]
% [tau1;tau2;...] ... sily/momenty v kloubech manipulatoru

jointCoords = varargin{1};
endEff_force_moment = varargin{2}; % w.r.t. Fb => feb, ueb
kinPar = varargin{3};
dynPar = varargin{4};

DHpar = kinPar.DHpar;
jointType = kinPar.jointType;
Qhome = kinPar.Qhome;

mass = dynPar.mass;
inertiaTensor = dynPar.inertiaTensor;
gravityCenter = dynPar.gravityCenter;
gravityVector = dynPar.gravityVector;
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

% gravitacni vektor v s.s. Fb
g0__0 = R0b'*gravityVector;

N = size(DHpar,1);

% transformacni matice
q = jointCoords(:,1);

% TT{i} = T_{i}^{i-1} => TT = {T10,T21,T32,T43,...}
% TT_0{i} = T_{i}^{0} => TT_0 = {T10,T20,T30,T40,...}
for i = 1:N
    if strcmp(jointType(i),'P')
        DHpar(i,1) = q(i);
        sigma(i) = 1;
    else
        DHpar(i,2) = q(i);
        sigma(i) = 0;   
    end
    sigma_hat(i) = 1 - sigma(i);
    
    TT{i} = DH(DHpar(i,:));
    
    if i > 1
        TT_0{i} = TT_0{i-1}*TT{i};
    else
        TT_0{i} = TT{i};
    end
end

% DOPREDNA REKURZE (rychlosti, zrychleni)
% pocatecni podminky
dO_vect = zeros(3,1); % dO_vect = [dO_0^0,dO_1^1,...];
omega_vect = zeros(3,1); % omega_vect = [omega_0^0,omega_1^1,...];
ddO_vect = zeros(3,1); % ddO_vect = [ddO_0^0,ddO_1^1,...];
domega_vect = zeros(3,1); % domega_vect = [domega_0^0,domega_1^1,...];

dq = jointCoords(:,2);
ddq = jointCoords(:,3);
for j = 1:N
    % rekurzivni vypocet rychlosti s.s. (translacni, uhlova) v s.s. akt.
    % ramene, r__a__b__c = r_{a,b}^c, R__a__b = R_a^b
    R__j_1__j = TT{j}(1:3,1:3)';
    r__j_1__j__j = R__j_1__j*TT{j}(1:3,4);

    omega_vect(:,j+1) = R__j_1__j*(omega_vect(:,j) + [0;0;sigma_hat(j)*dq(j)]);
    dO_vect(:,j+1) = R__j_1__j*(dO_vect(:,j) + [0;0;sigma(j)*dq(j)]) + cross(omega_vect(:,j+1),r__j_1__j__j);
    
    % rekurzivni vypocet zrychleni s.s. v s.s. akt. ramene
    domega_vect(:,j+1) = R__j_1__j*(domega_vect(:,j) + sigma_hat(j)*(dq(j)*[omega_vect(2,j);-omega_vect(1,j);0]+[0;0;ddq(j)]));
    ddO_vect(:,j+1) = R__j_1__j*(ddO_vect(:,j)+ dq(j)*sigma(j)*[omega_vect(2,j);-omega_vect(1,j);0] + sigma(j)*[0;0;ddq(j)]) + cross(domega_vect(:,j+1),r__j_1__j__j) + cross(omega_vect(:,j+1),cross(omega_vect(:,j+1),r__j_1__j__j)) + cross(omega_vect(:,j+1),sigma(j)*dq(j)*R__j_1__j(:,3));
    
    % zrychelni teziste ramen v s.s. akt. ramene: ddpC_vect = [ddpC_1^1,ddpC_2^2,...]
    r__j__Cj__j = gravityCenter(:,j);
    
    ddpC_vect(:,j) = ddO_vect(:,j+1) + cross(domega_vect(:,j+1),r__j__Cj__j) + cross(omega_vect(:,j+1),cross(omega_vect(:,j+1),r__j__Cj__j));
end

% ZPETNA REKURZE (sily, momenty)
% sila/moment f_i, u_i v s.s. akt. ramene: f_vect = [f_1^1,f_2^2,...,f_N^N,f_N+1^N],
% u_vect = [u_1^1,u_2^2,...,u_N^N,u_N+1^N]
% f_N+1^N, u_N+1^N ... sila, moment pusobici na konc. efektor w.r.t. s.s. F_N

% pro j = N+1: f_N+1^N = (R_N^0)^T*f_N+1^0, u_N+1^N = (R_N^0)^T*u_N+1^0,
% kde -f_N+1^0, -u_N+1^0 je sila/moment pusobici na konc. ef. vyjadreny
% vzhledem k s.s. zakladny F0

% f_i, u_i opdovidaji reakcnim silam/momentum v SimMechanicsu (nastaveni: Reactions mesured on: FOLLOWER, W.r.t. CS: LOCAL BODY)!

% kompenzace polohy zakladny a konc. eff
RN0 = TT_0{end}(1:3,1:3);
r__N__e__N = TeN(1:3,4);

% f_e^b, u_e^b -> f_N+1^0, u_N+1^0
f__Nplus1__0 = -R0b'*endEff_force_moment(1:3);
u__Nplus1__0 = -R0b'*(endEff_force_moment(4:6) + cross(R0b*RN0*r__N__e__N,endEff_force_moment(1:3)));

f_vect(:,N+1) = TT_0{N}(1:3,1:3)'*f__Nplus1__0; % f_N+1^N
u_vect(:,N+1) = TT_0{N}(1:3,1:3)'*u__Nplus1__0; % u_N+1^N 

for j = (N:-1:1)
    r__j__Cj__j = gravityCenter(:,j);
    R__j_1__j = TT{j}(1:3,1:3)';
    r__j_1__j__j = R__j_1__j*TT{j}(1:3,4);
    
    if j == N
        f_vect(:,j) = f_vect(:,j+1) + mass(j)*(ddpC_vect(:,j)-TT_0{j}(1:3,1:3)'*g0__0);
        u_vect(:,j) = u_vect(:,j+1) + cross(f_vect(:,j+1),r__j__Cj__j) - cross(f_vect(:,j),(r__j_1__j__j + r__j__Cj__j))+inertiaTensor(:,(3*(j-1)+1):(3*(j-1)+3))*domega_vect(:,j+1) + cross(omega_vect(:,j+1),inertiaTensor(:,(3*(j-1)+1):(3*(j-1)+3))*omega_vect(:,j+1));
    else
        f_vect(:,j) = TT{j+1}(1:3,1:3)*f_vect(:,j+1) + mass(j)*(ddpC_vect(:,j)-TT_0{j}(1:3,1:3)'*g0__0);
        u_vect(:,j) = TT{j+1}(1:3,1:3)*u_vect(:,j+1) + cross(TT{j+1}(1:3,1:3)*f_vect(:,j+1),r__j__Cj__j) - cross(f_vect(:,j),(r__j_1__j__j + r__j__Cj__j))+inertiaTensor(:,(3*(j-1)+1):(3*(j-1)+3))*domega_vect(:,j+1) + cross(omega_vect(:,j+1),inertiaTensor(:,(3*(j-1)+1):(3*(j-1)+3))*omega_vect(:,j+1));
    end
end

% PRUMETY REAKCNICH SIL (f_i) / MOMENTU (u_i) DO OS AKTUATORU => tau_i [N/Nm]
% tau_vect = [tau_1;tau_2;...]
for i = 1:size(f_vect,2)-1
    z__i_1__i = TT{i}(3,1:3)';
    tau_vect(i,1) = f_vect(:,i)'*z__i_1__i*sigma(i) + u_vect(:,i)'*z__i_1__i*sigma_hat(i);
end

varargout{1} = tau_vect;