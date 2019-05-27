% interpolace dat - externi feedrate
close all
clear
clc

% interpolovane body - translace
x = [1,2,3,-1,2,2+3*0.6]; %rand(5,1); %[1;2;-1;3];
y = [1,2,3,3,0,0-3*0.6]; %rand(5,1); %[1;2;-1;3];
z = [1,0,3,1,5,5+4*0.6];% rand(5,1); %[1;2;-1;3];
Pvect = [x;y;z];

% interpolovane body - orientace (kvaternion)
EA = [[0;0;0],[1;0.5;1],[1;1;0],[1;1;0;],[1;1;0],[0;0;0]];
% EA = [[1;1;0],[1;1;0],[1;1;0],[1;1;0;],[1;1;0]];

for i = 1:size(EA,2)
    Qvect(:,i) = rotMatrixAndAngularVel2Quaternion(EulerAnglesZYX2rotMatrixAndAngularVel(EA(:,i)));
end


% Radius vektor - parametr blendingu
radiusVect = [1;0.5;2;1];

% typ feedrate
feedType = 2;

data = linSegPolyBlendInterp_withOrient_parComp(Pvect,Qvect,radiusVect,feedType)



% profil rychlosti
S = data.S_knot(end)
v_max = 2;
a_max = 2;
b_max = 5;

% rozliseni generovani trajektoerie
res = 0.005;

[pos,vel,accel,jerk] = BAVS_genTrajectories(S,v_max,a_max,b_max,res);
time = pos(1,:);
S = pos(2,:);
V = vel(2,:);
A = accel(2,:);

% trajectory.signals.values(:,:,i) = [[O;Quat],[dO_dS;dQuat_dS],[d2O_dS2;d2Q_dS2]]
% feed.signals.values(:,:,i) = [u,du_dS,d2u_dS2]

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
    
    rotMangVel = Quaternion2rotMatrixAndAngularVel([Quat,dQuat_dt,d2Quat_dt2]);
    omega = rotMangVel(:,4);
    domega_dt = rotMangVel(:,5);
    
    trajectory.signals.values(:,:,i) = [[O,dO_dt,d2O_dt2];[Quat,dQuat_dt,d2Quat_dt2]];
    trajectory_omega.signals.values(:,:,i) = [omega,domega_dt];
    feed.signals.values(:,:,i) = [u,du_dt,d2u_dt2];
end


figure('Name','INTERPOLACE TRANSLACE')
subplot(3,1,1)
plot(trajectory.time,squeeze(trajectory.signals.values(1:3,1,:)))
xlabel('time')
ylabel('POS (xyz)')
subplot(3,1,2)
hold on
plot(trajectory.time,squeeze(trajectory.signals.values(1:3,2,:)))
plot(trajectory.time(1:end-1),[diff(squeeze(trajectory.signals.values(1,1,:)))'./diff(trajectory.time);diff(squeeze(trajectory.signals.values(2,1,:)))'./diff(trajectory.time);diff(squeeze(trajectory.signals.values(3,1,:)))'./diff(trajectory.time)],'--')
plot(trajectory.time,sqrt(squeeze(trajectory.signals.values(1,2,:)).^2 + squeeze(trajectory.signals.values(2,2,:)).^2 + squeeze(trajectory.signals.values(3,2,:)).^2),'LineWidth',2,'Color','k')
xlabel('time')
ylabel('VEL (xyz) + normVEL (xyz)')
subplot(3,1,3)
hold on
plot(trajectory.time,squeeze(trajectory.signals.values(1:3,3,:)))
plot(trajectory.time(1:end-1),[diff(squeeze(trajectory.signals.values(1,2,:)))'./diff(trajectory.time);diff(squeeze(trajectory.signals.values(2,2,:)))'./diff(trajectory.time);diff(squeeze(trajectory.signals.values(3,2,:)))'./diff(trajectory.time)],'--')
plot(trajectory.time,sqrt(squeeze(trajectory.signals.values(1,3,:)).^2 + squeeze(trajectory.signals.values(2,3,:)).^2 + squeeze(trajectory.signals.values(3,3,:)).^2),'LineWidth',2,'Color','k')
xlabel('time')
ylabel('ACCEL (xyz) + normACCEL (xyz)')

figure('Name','INTERPOLACE ORIENTACE')
subplot(3,1,1)
plot(trajectory.time,squeeze(trajectory.signals.values(4:7,1,:)))
xlabel('time')
ylabel('ORIENT (Quat)')
subplot(3,1,2)
hold on
plot(trajectory.time,squeeze(trajectory.signals.values(4:7,2,:)))
plot(trajectory.time(1:end-1),[diff(squeeze(trajectory.signals.values(4,1,:)))'./diff(trajectory.time);diff(squeeze(trajectory.signals.values(5,1,:)))'./diff(trajectory.time);diff(squeeze(trajectory.signals.values(6,1,:)))'./diff(trajectory.time);diff(squeeze(trajectory.signals.values(7,1,:)))'./diff(trajectory.time)],'--')
plot(trajectory_omega.time,sqrt(squeeze(trajectory_omega.signals.values(1,1,:)).^2 + squeeze(trajectory_omega.signals.values(2,1,:)).^2 + squeeze(trajectory_omega.signals.values(3,1,:)).^2),'LineWidth',2,'Color','k')
xlabel('time')
ylabel('VEL (Quat) + normOMEGA')
subplot(3,1,3)
hold on
plot(trajectory.time,squeeze(trajectory.signals.values(4:7,3,:)))
plot(trajectory.time(1:end-1),[diff(squeeze(trajectory.signals.values(4,2,:)))'./diff(trajectory.time);diff(squeeze(trajectory.signals.values(5,2,:)))'./diff(trajectory.time);diff(squeeze(trajectory.signals.values(6,2,:)))'./diff(trajectory.time);diff(squeeze(trajectory.signals.values(7,2,:)))'./diff(trajectory.time)],'--')
plot(trajectory_omega.time,sqrt(squeeze(trajectory_omega.signals.values(1,2,:)).^2 + squeeze(trajectory_omega.signals.values(2,2,:)).^2 + squeeze(trajectory_omega.signals.values(3,2,:)).^2),'LineWidth',2,'Color','k')
xlabel('time')
ylabel('ACCEL (Qaut) + normdOMEGA/dt')

figure('Name','FEEDRATE (pocitan z TRANSLACE!)')
hold on
plot(feed.time,squeeze(feed.signals.values(:,:,:)))
xlabel('time')
legend('u(t)','du(t)/dt','d^2u(t)/dt^2')


% animace translace + rotace
u_knot = data.u_knot;
S_knot = data.S_knot;
u_int = [];
for i = 1:size(Pvect,2)
    if i == 1
        u_int = [u_int,u_knot(1)];
    elseif i == size(Pvect,2)
        u_int = [u_int,u_knot(end)];
     else
        u_int = [u_int,u_knot(2*i-2)+(-u_knot(2*i-2)+u_knot(2*i-1))/2];
    end
end

for i = 1:length(S_knot)
    time_knot(i) = interp1(S,time,S_knot(i));
end

for i = 1:length(u_int)
    time_int(i) = interp1(u_knot,time_knot,u_int(i));
end

rotationAnimation(trajectory.time,squeeze(trajectory.signals.values(1:3,1,:)),squeeze(trajectory.signals.values(4:7,1,:)),time_int,Pvect,Qvect)

