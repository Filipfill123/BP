function [critFcnVal,critFcnVal_index] = critFcn(ksi,trajectory)

% optimalizovane parametry
% ksi = [L1;L2;L3];

% vyber reseni - mozno take optimalizovat... reseni IGM
sol = -1;

% endEff_force_moment
endEff_force_moment = [0;0;0;0;0;0];

% DH parametry
% DHpar = [d_1,theta_1,a_1,alpha_1; d_2,theta_2,a_2,alpha_2;... ]
DHpar = [0,0,ksi(1),0;0,0,ksi(2),0;0,0,ksi(3),0];

% typy kloubu
% jointType = [joint1, joint2,...] => kl. sour: [q1; q2; ...]
jointType = ['R','R','R'];

% domovska poloha - hodnoty kl. souradnic pro sestaveni modelu v SimMech
% Qhome = [q1h; q2h; ...]
Qhome = [0;0;0];

% kompenzace polohy koncoveho efektoru
% endEffComp = [Oen,Ren]
Ren = EulerAnglesXYZ2rotMatrix([0;0;0]);
Oen = [0;0;0];

endEffComp = [Oen,Ren];

% kompenzace polohy zakladny
% baseComp = [O0b,R0b]
R0b = EulerAnglesXYZ2rotMatrix([0;0;0]);
O0b = [0;0;0];

baseComp = [O0b,R0b];

% DYNAMICKE PARAMETRY MANIPULATORU + PARAMETRY POHONU
% Pozn.: bez uvazovani dynamickeho chovani rotoru motoru !!! (nekdy dodelat)

% model ramene robotu - plna tyc delky L, prumeru d, hustoty rho a s motorem o hmotnosti M
d = 0.03;
rho = 7800;

M1 = 1.5; % hmotnost 1. motoru
[m1,T1,I1] = simpleLinkDynPar(ksi(1),d,rho,M1);
M2 = 1; % hmotnost 2. motoru
[m2,T2,I2] = simpleLinkDynPar(ksi(2),d,rho,M2);
M3 = 0.5; % hmotnost 3. motoru
[m3,T3,I3] = simpleLinkDynPar(ksi(3),d,rho,M3);

% hmotnosti ramen Link1, Link2,...
% mass = [m1; m2; ...]
mass = [m1;m2;m3];

% moment setrvacnosti ramen Link1, Link2,... vyjadrena v s.s. daneho ramene
%   F1, F2, ...
%   => inertiaTensor = [I1, I2, ...] ... konstantni (nezavisi na poloze manip.!)
inertiaTensor = [I1,I2,I3];

% umisteni tezist ramen Link1, Link2, ... vyjadrena v s.s. daneho ramene
%   F1, F2, ...
%   => gravityCenter = [T1, T2, ...] ... konstantni (nezavisi na poloze
%   manip.!)
gravityCenter = [T1,T2,T3];


% gravitacni vektor v s.s. F0
gravityVector = [0;-9.81;0];

kinPar{1}.DHpar = DHpar;
kinPar{1}.jointType = jointType;
kinPar{1}.Qhome = Qhome;
kinPar{1}.endEffComp = endEffComp;
kinPar{1}.baseComp = baseComp;

dynPar{1}.mass = mass;
dynPar{1}.inertiaTensor = inertiaTensor;
dynPar{1}.gravityCenter = gravityCenter;
dynPar{1}.gravityVector = gravityVector;

% silove momenty aktuatoru podel trajektorie
for i = 1:length(trajectory.time)
    jointCoords = inverseKinematics(trajectory.signals.values(:,:,i),kinPar{1},sol);
    tau(:,i) = inverseDynamicModel(jointCoords,endEff_force_moment,kinPar{1},dynPar{1});
    tau_norm(i) = norm(tau(:,i));
end

% mira optimality (max. moment, prum. moment,...)
[critFcnVal,critFcnVal_index] = max(tau_norm); % maximalni norma momentu podel trajektorie


