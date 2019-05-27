close all hidden
clear
clc
addpath('genTraj')

% **************************************
% ZADANI PARAMETRU - ZACATEK
% **************************************

% KINEMATICKE PARAMETRY MANIPULATORU

% optimalizovane parametry
L1 = 1;
L2 = 1;
L3 = 1;

ksi = [L1;L2;L3];

% vyber reseni
sol = -1;

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
[m1,T1,I1] = simpleLinkDynPar(L1,d,rho,M1);
M2 = 1; % hmotnost 2. motoru
[m2,T2,I2] = simpleLinkDynPar(L2,d,rho,M2);
M3 = 0.5; % hmotnost 3. motoru
[m3,T3,I3] = simpleLinkDynPar(L3,d,rho,M3);

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


% GENERATOR TRAJEKTORIE - linearni segmenty s polynomialnim napojeni

% interpolovane body
x_vect = [1,1.3,1.3,1,1];
y_vect = [1,1,1.3,1.3,1];
phi_vect = [0,pi/8,pi/4,pi/8,0];

% koincidencni body v translaci Z = 0
Pvect = [x_vect;y_vect;0*x_vect];

% koincidenci body v orientaci, jen rotace kolem Z
for i = 1:size(phi_vect,2) % rotace v ose Z o uhel phi_vect vyjadrena v Kvaternionu - pro interpolator
    Qvect(:,i) = rotMatrixAndAngularVel2Quaternion(EulerAnglesZYX2rotMatrixAndAngularVel([phi_vect(i);0;0]));
end

% Radius vektor - parametr blendingu
% radiusVect = 1*ones(1,size(Pvect,2)-2);
radiusVect = [0.03;0.03;0.03];

% typ feedrate
feedType = 1;

% interpolator
data = linSegPolyBlendInterp_withOrient_parComp(Pvect,Qvect,radiusVect,feedType);

% profil rychlosti
S = data.S_knot(end); % celkove ujeta draha
v_max = 0.2; % max. rychlost
a_max = 1.5; % max. zrychleni

% rozliseni generovani trajektoerie
res = 0.05;

[pos,vel,accel] = PVA_genTrajectories(S,v_max,a_max,res);
time = pos(1,:);
S = pos(2,:);
V = vel(2,:);
A = accel(2,:);

trajectory.time = time;
trajectory_omega.time = time;
feed.time = time;
for i = 1:length(S)
    [trans,orient,FEED] = linSegPolyBlendInterp_withOrient_eval(S(i),data);
    O = trans(:,1);
    dO_dS = trans(:,2);
    d2O_dS2 = trans(:,3);
    
    Quat = orient(:,1);
    dQuat_dS = orient(:,2);
    d2Quat_dS2 = orient(:,3);
    
    u = FEED(:,1);
    du_dS = FEED(:,2);
    d2u_dS2 = FEED(:,3);
    
    % zavislosti na case
    dO_dt = dO_dS*V(i);
    d2O_dt2 = d2O_dS2*V(i)^2 + dO_dS*A(i);
    dQuat_dt = dQuat_dS*V(i);
    d2Quat_dt2 = d2Quat_dS2*V(i)^2 + dQuat_dS*A(i);
    du_dt = du_dS*V(i);
    d2u_dt2 = d2u_dS2*V(i)^2 + du_dS*A(i);
    
%     rotMangVel = Quaternion2rotMatrixAndAngularVel([Quat,dQuat_dt,d2Quat_dt2]);
%     omega = rotMangVel(:,4);
%     domega_dt = rotMangVel(:,5);
    EA = rotMatrixAndAngularVel2EulerAnglesZYX(Quaternion2rotMatrixAndAngularVel([Quat,dQuat_dt,d2Quat_dt2]));
    trajectory.signals.values(:,:,i) = [[O(1:2),dO_dt(1:2),d2O_dt2(1:2)];EA(1,:)];
%     trajectory_omega.signals.values(:,:,i) = [omega,domega_dt];
    feed.signals.values(:,:,i) = [u,du_dt,d2u_dt2];
end


figure('Name','INTERPOLACE TRANSLACE')
subplot(3,1,1)
plot(trajectory.time,squeeze(trajectory.signals.values(1:2,1,:)))
xlabel('time')
ylabel('POS (xyz)')
subplot(3,1,2)
hold on
plot(trajectory.time,squeeze(trajectory.signals.values(1:2,2,:)))
plot(trajectory.time(1:end-1),[diff(squeeze(trajectory.signals.values(1,1,:)))'./diff(trajectory.time);diff(squeeze(trajectory.signals.values(2,1,:)))'./diff(trajectory.time)],'--')
plot(trajectory.time,sqrt(squeeze(trajectory.signals.values(1,2,:)).^2 + squeeze(trajectory.signals.values(2,2,:)).^2),'LineWidth',2,'Color','k')
xlabel('time')
ylabel('VEL (xyz) + normVEL (xyz)')
subplot(3,1,3)
hold on
plot(trajectory.time,squeeze(trajectory.signals.values(1:2,3,:)))
plot(trajectory.time(1:end-1),[diff(squeeze(trajectory.signals.values(1,2,:)))'./diff(trajectory.time);diff(squeeze(trajectory.signals.values(2,2,:)))'./diff(trajectory.time)],'--')
plot(trajectory.time,sqrt(squeeze(trajectory.signals.values(1,3,:)).^2 + squeeze(trajectory.signals.values(2,3,:)).^2),'LineWidth',2,'Color','k')
xlabel('time')
ylabel('ACCEL (xyz) + normACCEL (xyz)')

figure('Name','INTERPOLACE ORIENTACE')
subplot(3,1,1)
plot(trajectory.time,squeeze(trajectory.signals.values(3,1,:)))
xlabel('time')
ylabel('ORIENT (phi)')
subplot(3,1,2)
hold on
plot(trajectory.time,squeeze(trajectory.signals.values(3,2,:)))
plot(trajectory.time(1:end-1),diff(squeeze(trajectory.signals.values(3,1,:)))'./diff(trajectory.time),'--')
xlabel('time')
ylabel('VEL (phi)')
subplot(3,1,3)
hold on
plot(trajectory.time,squeeze(trajectory.signals.values(3,3,:)))
plot(trajectory.time(1:end-1),diff(squeeze(trajectory.signals.values(3,2,:)))'./diff(trajectory.time),'--')
xlabel('time')
ylabel('ACCEL (phi)')

figure('Name','FEEDRATE (pocitan z TRANSLACE!)')
hold on
plot(feed.time,squeeze(feed.signals.values(:,:,:)))
xlabel('time')
legend('u(t)','du(t)/dt','d^2u(t)/dt^2')

% uprava - neni souradnice z
trajectory.signals.dimensions = [3,3];
trajectory.signals.values = trajectory.signals.values([1:2,3],:,:);

figure('Name','Interpolace XYZ')
hold on
plot(squeeze(trajectory.signals.values(1,1,:)),squeeze(trajectory.signals.values(2,1,:)),'LineWidth',3)
plot(x_vect,y_vect,'Marker','.','MarkerSize',20)
xlabel('x')
ylabel('y')


% **************************************
% ZADANI PARAMETRU - KONEC
% **************************************

for i = 1:size(Qhome)
    if strcmp(jointType(i),'P')
        DHpar(i,1) = Qhome(i);
    else
        DHpar(i,2) = Qhome(i);
    end
end

kinPar{1}.DHpar = DHpar;
kinPar{1}.jointType = jointType;
kinPar{1}.Qhome = Qhome;
kinPar{1}.endEffComp = endEffComp;
kinPar{1}.baseComp = baseComp;

dynPar{1}.mass = mass;
dynPar{1}.inertiaTensor = inertiaTensor;
dynPar{1}.gravityCenter = gravityCenter;
dynPar{1}.gravityVector = gravityVector;

Links{1} = robotSetup(kinPar{1},dynPar{1});


