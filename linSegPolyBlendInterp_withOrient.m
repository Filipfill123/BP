function [trajectory,arcLength,feed,trajectory_omega,trajectory_u,u_knot,time_knot,S_knot] = linSegPolyBlendInterp_withOrient(Pvect,Qvect,radiusVect,v_max,a_max,b_max,res,feedType)
% Pvect = [P1,P2,...], Qvect = [Q1,Q2,...] ... interpolovane body translace, orientace (kvaternion) (feedrate se pocita pouze z translace)
% radiusVect ... vektor radiusu v translaci, kdy zacina blending na linearnich segmentech (vzd. okamziku blendovani od vrcholu Pi)
%   Interpolace orientace v lin. segmentech = SLERP (od Ovect do Ivect - pocatek konec blendovani v translaci) + prima interpolace prvku kvaternionu (ne zcela dobre) v blending segmentech kubickym SPLINEM (od Ivect do Ovect)
% v_max, a_max, b_max ... omezeni na rychlost, zrychleni, jerk podel trajektorie (pro feedrate) - model BAVS, pokud b_max = nan - jen model AVS
% res ... rozliseni v [s] casoveho okamziku generovani trajektorie
% feedType  = 1 (aproximativni vypocet feedrate = fluktuace tecne rychlosti, zpomaleni v beldning segmentech - jednoduchy ALG.)
%           = 2 (dodrzeni predepsaneho profilu rychlosti BAVS podel trajektorie, numer. ALG)
% 
% trajectory.signals.values(:,:,i) = [[O(t);Quat(t)],[dO(t)_dt;dQuat(t)_dt],[d2O(t)_dt2;d2Q(t)_dt2]]
% feed.signals.values(:,:,i) = [u(t),du(t)_dt,d2u(t)_dt2]
% arcLength.signals.values(:,:,i) = [s(t),v(t),a(t)]
% trajectory_omega.signals.values(:,:,i) = [omega(t),domega(t)_dt]

numOfPoints = size(Pvect,2);

% POZOR: osetreni stejnych koincidencnich bodu (spravne by to bylo osetrit dale v algoritmu! ZDE: pokud jsou dva po sobe jdouci koincidencni body shodne, jeden se perturbuje malou hodnotou)
for i = 1:size(Pvect,2)-1
    % translace
    if norm(Pvect(:,i)-Pvect(:,i+1)) < 1e-8
        disp('******** IDENTICKE ZADANE TRANSLACE !!! **********')
        disp(['Pvect(:,',num2str(i),') = Pvect(:,',num2str(i+1),') = [',num2str(Pvect(:,i)'),']'])
        % mala perturbace
        Pvect(:,i+1) = Pvect(:,i+1) + 1e-5*rand(3,1);
        disp([' => (perturbace) Pvect(:,',num2str(i+1),') = [',num2str(Pvect(:,i+1)'),']'])
        disp('******************')
    end
    % orientace
    dQuat = quatMult(Qvect(:,i),quatInv(Qvect(:,i+1)));
    if abs(abs(dQuat(1))-1) < 1e-8
        disp('******** IDENTICKE ZADANE ORIENTACE !!! **********')
        disp(['Qvect(:,',num2str(i),') = Qvect(:,',num2str(i+1),') = [',num2str(Qvect(:,i)'),']'])
        % mala perturbace
        Qvect(:,i+1) = Qvect(:,i+1) + 1e-5*rand(4,1);
        Qvect(:,i+1) = Qvect(:,i+1)/norm(Qvect(:,i+1));
        disp([' => (perturbace) Qvect(:,',num2str(i+1),') = [',num2str(Qvect(:,i+1)'),']'])
        disp('******************')
    end
    
end

% osetreni spatne nastavenych radiusu
for i = 2:length(radiusVect)-1
    if norm(Pvect(:,i+1)-Pvect(:,i)) <= radiusVect(i-1) + radiusVect(i)
        error(['Radius ',num2str(i-1),' a/nebo ',num2str(i),' je prilis velky.'])
    end
end
if norm(Pvect(:,2)-Pvect(:,1)) <= radiusVect(1);
    error(['Prvni radius je prilis velky.'])
end
if norm(Pvect(:,end)-Pvect(:,end-1)) <= radiusVect(end);
    error(['Posledni radius je prilis velky.'])
end

% pocatky a konce blendovani
% translace
for i = 2:numOfPoints-1
    % translace
    Ovect(:,i-1) = (Pvect(:,i+1)-Pvect(:,i))/norm(Pvect(:,i+1)-Pvect(:,i))*radiusVect(i-1) + Pvect(:,i);
    Ivect(:,i-1) = -(Pvect(:,i)-Pvect(:,i-1))/norm(Pvect(:,i)-Pvect(:,i-1))*radiusVect(i-1) + Pvect(:,i);
end

% orientace (parametr theta pro zacatek a konec blendovani (pro Ivect, Ovect) vyplyva z SLERP inter. mezi zadanymi kvaterniony v koincidencnich bodech)
% interpolace v lin. segmentech = SLERP (od Ovect do Ivect), prima interpolace prvku kvaternionu (ne zcela dobre) v blending segmentech kubickym SPLINEM (od Ivect do Ovect)
for i = 1:numOfPoints-1
    deltaTheta(:,i) = 2*acos(Qvect(:,i)'*Qvect(:,i+1));
    if (abs(deltaTheta(:,i)) > pi) % prejezd "kratsi" stranou
        Qvect(:,i+1) = -Qvect(:,i+1);
        deltaTheta(:,i) = 2*acos(Qvect(:,i)'*Qvect(:,i+1));
    end
end

for i = 1:numOfPoints-2
    deltaThetaIvect(:,i) = (1-norm(Pvect(:,i+1)-Ivect(:,i))/norm(Pvect(:,i+1)-Pvect(:,i)))*deltaTheta(:,i); % delka o kterou ujede theta na linearnim segmentu
    deltaThetaOvect(:,i) = (norm(Pvect(:,i+1)-Ovect(:,i))/norm(Pvect(:,i+2)-Pvect(:,i+1)))*deltaTheta(:,i+1);
end

% delky segmentu - APROXIMACE
arcLength_polySeg_approx = polySegLength_approx(Pvect,Ovect,Ivect);
arcLength_linSeg = linSegLength(Pvect,Ovect,Ivect);

% APROXIMACE: ujeta draha + uzlovy vektor u = <0,1>
S_knot = [0];
for i = 1:numOfPoints-2
    S_knot = [S_knot, S_knot(end) + arcLength_linSeg(i)];
    S_knot = [S_knot,S_knot(end) + arcLength_polySeg_approx(i)];
end
S_knot_approx = [S_knot,S_knot(end) + arcLength_linSeg(end)];
        
S_max_approx = S_knot_approx(end);
u_knot = S_knot_approx/S_max_approx;

[linSegPar_trans,linSegPar_orient] = linSegParam(u_knot,Pvect,Ovect,Ivect,deltaThetaOvect,deltaThetaIvect,deltaTheta,Qvect);
[polySegPar_trans,polySegPar_orient] = polySegParam(u_knot,linSegPar_trans,linSegPar_orient,Qvect,deltaTheta);

% vypocet feedrate (dle TRANSLACE!!!!!) + generovani trajektorie
switch feedType
    case 1 % aproximace
        S_knot = S_knot_approx;
        
        % generator profilu pohybu po trajektorii
        if isnan(b_max)
            [pos,vel,accel] = PVA_genTrajectories(S_max_approx,v_max,a_max,res);
            time = pos(1,:);
            S = pos(2,:);
            V = vel(2,:);
            A = accel(2,:);
        else
            [pos,vel,accel,jerk] = BAVS_genTrajectories(S_max_approx,v_max,a_max,b_max,res);
            time = pos(1,:);
            S = pos(2,:);
            V = vel(2,:);
            A = accel(2,:);
            B = jerk(2,:);
        end
        
        for i = 1:length(time)
            u(1,i) = S(i)/S_max_approx;
            du_dt(1,i) = V(i)/S_max_approx;
            d2u_dt2(1,i) = A(i)/S_max_approx;
        end
        
    case 2 % presne
        
        % skutecna delka polynomialnich segmentu
        S_knot = [0];
        arcLength_polySeg = polySegLength(u_knot,polySegPar_trans);
        for i = 1:numOfPoints-2
            S_knot = [S_knot, S_knot(end) + arcLength_linSeg(i)];
            S_knot = [S_knot,S_knot(end) + arcLength_polySeg(i)];
        end
        S_knot = [S_knot,S_knot(end) + arcLength_linSeg(end)];
        S_max = S_knot(end);
        
        % generator profilu pohybu po trajektorii
        if isnan(b_max)
            [pos,vel,accel] = PVA_genTrajectories(S_max,v_max,a_max,res);
            time = pos(1,:);
            S = pos(2,:);
            V = vel(2,:);
            A = accel(2,:);
        else
            [pos,vel,accel,jerk] = BAVS_genTrajectories(S_max,v_max,a_max,b_max,res);
            time = pos(1,:);
            S = pos(2,:);
            V = vel(2,:);
            A = accel(2,:);
            B = jerk(2,:);
        end

        for i = 1:length(time)
            FEED = feedrate([S(i);V(i);A(i)],u_knot,S_knot,polySegPar_trans);
            u(1,i) = FEED(1);
            du_dt(1,i) = FEED(2);
            d2u_dt2(1,i) = FEED(3);
        end
    otherwise
        error('Spatny feedType')
end

% time_knot
for i = 1:length(S_knot)
    if i == length(S_knot)
        time_knot(i) = time(end);
    else
        time_knot(i) = interp1(S,time,S_knot(i));
    end
end

trajectory.time = time;
trajectory.signals.dimensions = [7,3]; % [x;y;z;Quat1;Quat2;Quat3;Quat4]

trajectory_u.time = time;
trajectory_u.signals.dimensions = [7,3]; % [x;y;z;Quat1;Quat2;Quat3;Quat4]

feed.time = time;
feed.signals.dimensions = [1,3];

trajectory_omega.time = time;
trajectory_omega.signals.dimensions = [3,2];

arcLength.time = time;
arcLength.signals.dimensions = [1,3];

for i = 1:length(u)
    [trans,orient] = linSegPolyBlendEval_withOrient(u(i),u_knot,linSegPar_trans,linSegPar_orient,polySegPar_trans,polySegPar_orient,deltaTheta);
   
    p(:,i) = trans(:,1);
    dp_du(:,i) = trans(:,2);
    d2p_du2(:,i) = trans(:,3);
    
    q(:,i) = orient(:,1);
    dq_du(:,i) = orient(:,2);
    d2q_du2(:,i) = orient(:,3);
    
    % casove zavislosti
    dp_dt(:,i) = dp_du(:,i)*du_dt(:,i);
    d2p_dt2(:,i) = d2p_du2(:,i)*du_dt(:,i)^2 + dp_du(:,i)*d2u_dt2(:,i);
    
    dq_dt(:,i) = dq_du(:,i)*du_dt(:,i);
    d2q_dt2(:,i) = d2q_du2(:,i)*du_dt(:,i)^2 + dq_du(:,i)*d2u_dt2(:,i);
    
    % translace: tecna rychlost
    dpNorm_dt(:,i) = norm(dp_dt(1:3,i));
    
    % rotace: uhlova rychlost a zrychleni
    rotMatVelAccel = Quaternion2rotMatrixAndAngularVel([q(:,i),dq_dt(:,i),d2q_dt2(:,i)]);
    R{i} = rotMatVelAccel(:,1:3);
    omega(:,i) = rotMatVelAccel(:,4);
    domega_dt(:,i) = rotMatVelAccel(:,5);
    
    trajectory.signals.values(:,:,i) = [[p(:,i);q(:,i)],[dp_dt(:,i);dq_dt(:,i)],[d2p_dt2(:,i);d2q_dt2(:,i)]];
    trajectory_u.signals.values(:,:,i) = [[p(:,i);q(:,i)],[dp_du(:,i);dq_du(:,i)],[d2p_du2(:,i);d2q_du2(:,i)]];
    feed.signals.values(:,:,i) = [u(:,i),du_dt(:,i),d2u_dt2(:,i)];
    arcLength.signals.values(:,:,i) = [S(i),V(i),A(i)];
    trajectory_omega.signals.values(:,:,i) = [omega(:,i), domega_dt(:,i)];
end

%{
% ************************
% VYKRESLENI
% ************************
col = ['bgrcmy'];

figure('Name','Feedrate - BAVS')
hold on
plot(time,S,col(1))
plot(time,V,col(2))
plot(time,A,col(3))
% plot(time,B,col(4))
xlabel('time')
legend('S','V','A')

figure('Name','Interpolace v XYZ: TRANSLACE')
hold on
plot3(p(1,:),p(2,:),p(3,:),'LineWidth',2,'Color','k')
plot3(Pvect(1,:),Pvect(2,:),Pvect(3,:),'Marker','o','MarkerSize',5)
plot3(Ivect(1,:),Ivect(2,:),Ivect(3,:),'Marker','.','MarkerSize',15,'LineStyle','none','Color','m')
plot3(Ovect(1,:),Ovect(2,:),Ovect(3,:),'Marker','.','MarkerSize',15,'LineStyle','none','Color','m')
axis equal
xlabel('X')
xlabel('Y')
xlabel('Z')

figure('Name','Interpolace v par u: TRANSLACE')
subplot(3,1,1)
hold on
for i = 1:size(p,1)
    plot(u,p(i,:),'Color',col(i))
end
xlabel('u')
ylabel('p(u)')
subplot(3,1,2)
hold on
for i = 1:size(dp_du,1)
    plot(u,dp_du(i,:),'Color',col(i))
    plot(u(1:end-1),[diff(p(1,:))./diff(u);diff(p(2,:))./diff(u);diff(p(3,:))./diff(u)],'--')
end
xlabel('u')
ylabel('dp(u)/du')
subplot(3,1,3)
hold on
for i = 1:size(d2p_du2,1)
    plot(u,d2p_du2(i,:),'Color',col(i))
    plot(u(1:end-1),[diff(dp_du(1,:))./diff(u);diff(dp_du(2,:))./diff(u);diff(dp_du(3,:))./diff(u)],'--')
end
xlabel('u')
ylabel('d^2p(u)/du^2')

figure('Name','Interpolace v case: TRANSLACE')
subplot(4,1,1)
hold on
for i = 1:size(p,1)
    plot(time,p(i,:),'Color',col(i))
end
xlabel('t')
ylabel('p(t)')
subplot(4,1,2)
hold on
for i = 1:size(dp_dt,1)
    plot(time,dp_dt(i,:),'Color',col(i))
    plot(time(1:end-1),diff(p(i,:))./diff(time),'Color',col(i),'LineStyle','--')
end
xlabel('t')
ylabel('dp(t)/dt')
subplot(4,1,3)
hold on
for i = 1:size(d2p_dt2,1)
    plot(time,d2p_dt2(i,:),'Color',col(i))
   plot(time(1:end-1),diff(dp_dt(i,:))./diff(time),'Color',col(i),'LineStyle','--')
end
xlabel('t')
ylabel('d^2p(t)/dt^2')

subplot(4,1,4)
plot(time,dpNorm_dt)
xlabel('t')
ylabel('norm(d^p(t)/dt)')

figure('Name','Interpolace v par u: ORIENTACE (kvaterniony)')
subplot(3,1,1)
plot(u,q)
xlabel('u')
ylabel('Q')
subplot(3,1,2)
hold on
plot(u,dq_du)
plot(u(1:end-1),[diff(q(1,:))./diff(u);diff(q(2,:))./diff(u);diff(q(3,:))./diff(u);diff(q(4,:))./diff(u)],'--')
xlabel('u')
ylabel('dQ/du')

subplot(3,1,3)
hold on
plot(u,d2q_du2)
plot(u(1:end-1),[diff(dq_du(1,:))./diff(u);diff(dq_du(2,:))./diff(u);diff(dq_du(3,:))./diff(u);diff(dq_du(4,:))./diff(u)],'--')
xlabel('u')
ylabel('d^2Q/du^2')

figure('Name','Interpolace v case: ORIENTACE (kvaterniony)')
subplot(3,1,1)
plot(time,q)
xlabel('time')
ylabel('Q')
subplot(3,1,2)
hold on
plot(time,dq_dt)
plot(time(1:end-1),[diff(q(1,:))./diff(time);diff(q(2,:))./diff(time);diff(q(3,:))./diff(time);diff(q(4,:))./diff(time)],'--')
xlabel('time')
ylabel('dQ/dt')
subplot(3,1,3)
hold on
plot(time,d2q_dt2)
plot(time(1:end-1),[diff(dq_dt(1,:))./diff(time);diff(dq_dt(2,:))./diff(time);diff(dq_dt(3,:))./diff(time);diff(dq_dt(4,:))./diff(time)],'--')
xlabel('time')
ylabel('d^2Q/dt^2')

figure('Name','Interpolace v case: ORIENTACE: Uhlova rychlost, zrychleni')
subplot(2,1,1)
plot(time,omega)
xlabel('time')
ylabel('omega')
subplot(2,1,2)
hold on
plot(time,domega_dt)
plot(time(1:end-1),[diff(omega(1,:))./diff(time);diff(omega(2,:))./diff(time);diff(omega(3,:))./diff(time)],'--')
xlabel('time')
ylabel('domega/dt')

figure('Name','Feedrate')
hold on
plot(time,u,'Color',col(1))
plot(time,du_dt,'Color',col(2))
plot(time,d2u_dt2,'Color',col(3))
xlabel('time')
legend('u(t)','du(t)/dt','du^2(t)/dt^2')
plot(time(1:end-1),diff(u)./diff(time),'Color',col(2),'LineStyle','--')
plot(time(1:end-1),diff(du_dt)./diff(time),'Color',col(3),'LineStyle','--')

%}