function [trajectory,feed] = linSegPolyBlendInterp(Pvect,radiusVect,v_max,a_max,b_max,res,feedType)
% Pvect = [P1,P2,...] interpolovane body (feedrate se pocita pouze z prvnich trech souradnic - typicky XYZ polohy)
% radiusVect ... vektor radiusu, kdy zacina blending na linearnich segmentech (vzd. okamziku blendovani od vrcholu Pi)
% v_max, a_max, b_max ... omezeni na rychlost, zrychleni, jerk podel trajektorie (pro feedrate) - model BAVS
% res ... rozliseni v [s] casoveho okamziku generovani trajektorie
% feedType  = 1 (aproximativni vypocet feedrate = fluktuace tecne rychlosti, zpomaleni v beldning segmentech - jednoduchy ALG.)
%           = 2 (dodrzeni predepsaneho profilu rychlosti BAVS podel trajektorie, numer. ALG)

numOfPoints = size(Pvect,2);

% osetreni spatne nastavenych radiusu
for i = 1:length(radiusVect)-1
    if norm(Pvect(:,i+2)-Pvect(:,i+1)) <= radiusVect(i) + radiusVect(i+1)
        error(['Radius ',num2str(i),' a/nebo ',num2str(i+1),' je prilis velky.'])
    end
end
if norm(Pvect(:,2)-Pvect(:,1)) <= radiusVect(1);
    error(['Prvni radius je prilis velky.'])
end
if norm(Pvect(:,end)-Pvect(:,end-1)) <= radiusVect(end);
    error(['Posledni radius je prilis velky.'])
end

% pocatky a konce blendovani
for i = 2:numOfPoints-1
    Ovect(:,i-1) = (Pvect(:,i+1)-Pvect(:,i))/norm(Pvect(:,i+1)-Pvect(:,i))*radiusVect(i-1) + Pvect(:,i);
    Ivect(:,i-1) = -(Pvect(:,i)-Pvect(:,i-1))/norm(Pvect(:,i)-Pvect(:,i-1))*radiusVect(i-1) + Pvect(:,i);
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

linSegPar = linSegParam(u_knot,Pvect,Ovect,Ivect);
polySegPar = polySegParam(u_knot,linSegPar);


% vypocet feedrate + generovani trajektorie
switch feedType
    case 1 % aproximace
        
        % generator profilu pohybu po trajektorii
        [pos,vel,accel,jerk] = BAVS_genTrajectories(S_max_approx,v_max,a_max,b_max,res);
        time = pos(1,:);
        S = pos(2,:);
        V = vel(2,:);
        A = accel(2,:);
        B = jerk(2,:);
        for i = 1:length(time)
            u(1,i) = S(i)/S_max_approx;
            du_dt(1,i) = V(i)/S_max_approx;
            d2u_dt2(1,i) = A(i)/S_max_approx;
        end
        
    case 2 % presne
        
        % skutecna delka polynomialnich segmentu
        S_knot = [0];
        arcLength_polySeg = polySegLength(u_knot,polySegPar);
        for i = 1:numOfPoints-2
            S_knot = [S_knot, S_knot(end) + arcLength_linSeg(i)];
            S_knot = [S_knot,S_knot(end) + arcLength_polySeg(i)];
        end
        S_knot = [S_knot,S_knot(end) + arcLength_linSeg(end)];
        S_max = S_knot(end);
        
        % generator profilu pohybu po trajektorii
        [pos,vel,accel,jerk] = BAVS_genTrajectories(S_max,v_max,a_max,b_max,res);
        time = pos(1,:);
        S = pos(2,:);
        V = vel(2,:);
        A = accel(2,:);
        B = jerk(2,:);

        for i = 1:length(time)
            feed = feedrate([S(i);V(i);A(i)],u_knot,S_knot,polySegPar);
            u(1,i) = feed(1);
            du_dt(1,i) = feed(2);
            d2u_dt2(1,i) = feed(3);
        end
    otherwise
        error('Spatny feedType')
end
        

trajectory.time = time;
trejectory.signals.dimensions = [size(Pvect,1),3];

feed.time = time;
feed.signals.dimensions = [1,3];

for i = 1:length(u)
    [p_int(:,i),dp_int_du(:,i),d2p_int_du2(:,i)] = linSegPolyBlendEval(u(i),u_knot,linSegPar,polySegPar);
    
    % casove zavislslosti
    dp_int_dt(:,i) = dp_int_du(:,i)*du_dt(:,i);
    d2p_int_dt2(:,i) = d2p_int_du2(:,i)*du_dt(:,i)^2 + dp_int_du(:,i)*d2u_dt2(:,i);
    
    % tecna rychlost
    dp_intNorm_dt(:,i) = norm(dp_int_dt(1:3,i));
    
    trajectory.signals.values(:,:,i) = [p_int(:,i),dp_int_dt(:,i),d2p_int_dt2(:,i)];
    feed.signals.values(:,:,i) = [u(:,i),du_dt(:,i),d2u_dt2(:,i)];
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
plot(time,B,col(4))
xlabel('time')
legend('S','V','A','B')

figure('Name','Interpolace v par u')
subplot(3,1,1)
hold on
for i = 1:size(p_int,1)
    plot(u,p_int(i,:),'Color',col(i))
end
xlabel('u')
ylabel('p(u)')
subplot(3,1,2)
hold on
for i = 1:size(dp_int_du,1)
    plot(u,dp_int_du(i,:),'Color',col(i))
end
xlabel('u')
ylabel('dp(u)/du')
subplot(3,1,3)
hold on
for i = 1:size(d2p_int_du2,1)
    plot(u,d2p_int_du2(i,:),'Color',col(i))
end
xlabel('u')
ylabel('d^2p(u)/du^2')

figure('Name','Interpolace v case')
subplot(4,1,1)
hold on
for i = 1:size(p_int,1)
    plot(time,p_int(i,:),'Color',col(i))
end
xlabel('t')
ylabel('p(t)')
subplot(4,1,2)
hold on
for i = 1:size(dp_int_dt,1)
    plot(time,dp_int_dt(i,:),'Color',col(i))
    plot(time(1:end-1),diff(p_int(i,:))./diff(time),'Color',col(i),'LineStyle','--')
end
xlabel('t')
ylabel('dp(t)/dt')
subplot(4,1,3)
hold on
for i = 1:size(d2p_int_dt2,1)
    plot(time,d2p_int_dt2(i,:),'Color',col(i))
   plot(time(1:end-1),diff(dp_int_dt(i,:))./diff(time),'Color',col(i),'LineStyle','--')
end
xlabel('t')
ylabel('d^2p(t)/dt^2')

subplot(4,1,4)
plot(time,dp_intNorm_dt)
xlabel('t')
ylabel('norm(d^p(t)/dt)')

figure('Name','Interpolace v XYZ')
hold on
plot3(p_int(1,:),p_int(2,:),p_int(3,:),'LineWidth',2,'Color','k')
plot3(Pvect(1,:),Pvect(2,:),Pvect(3,:),'Marker','o','MarkerSize',5)
plot3(Ivect(1,:),Ivect(2,:),Ivect(3,:),'Marker','.','MarkerSize',15,'LineStyle','none','Color','m')
plot3(Ovect(1,:),Ovect(2,:),Ovect(3,:),'Marker','.','MarkerSize',15,'LineStyle','none','Color','m')
axis equal
xlabel('X')
xlabel('Y')
xlabel('Z')


figure('Name','Feedrate')
hold on
plot(time,u,'Color',col(1))
plot(time,du_dt,'Color',col(2))
plot(time,d2u_dt2,'Color',col(3))
xlabel('time')
legend('u(t)','du(t)/dt','du^2(t)/dt^2')
plot(time(1:end-1),diff(u)./diff(time),'Color',col(2),'LineStyle','--')
plot(time(1:end-1),diff(du_dt)./diff(time),'Color',col(3),'LineStyle','--')

polySegPar{:}
%}
