function [pos,vel,accel,jerk] = BAVS_genTrajectories(S,maxVelocity,maxAcceleration,maxJerk,res)

% [pos,vel,accel] = PVA_genTrajectories(S,maxVelocity,maxAcceleration,res)
%
% generator primkoveho pohybu po draze (S) s omezenim na rychlost
% (maxVelocity) a zrychleni (maxAcceleration)
% res =  interval mezi casovymmi vzorky (priblizne)

V = abs(maxVelocity);
A = abs(maxAcceleration);
B = abs(maxJerk);

T = [0,BAVS(S,maxVelocity,maxAcceleration,maxJerk)];
jerk = [B,0,-B,0,-B,0,B];

s_T = 0;
v_T = 0;
a_T = 0;

time = [];
s = [];
v = [];
a = [];
b = [];
for i = 1:length(T)-1
    int{i} = linspace(T(i),T(i+1),round((T(i+1)-T(i))/res));
    
    if length(int{i}) > 0
        time_int = int{i} - int{i}(1);
        
        b = [b(1:end-1),jerk(i)*ones(size(time_int))];
        a = [a(1:end-1),a_T(i) + jerk(i)*time_int];
        v = [v(1:end-1),v_T(i) + a_T(i)*time_int + 1/2*jerk(i)*time_int.^2];
        s = [s(1:end-1),s_T(i) + v_T(i)*time_int + 1/2*a_T(i)*time_int.^2 + 1/6*jerk(i)*time_int.^3];
        
        time = [time(1:end-1),int{i}];
    end
    a_T(i+1) = a(end);
    v_T(i+1) = v(end);
    s_T(i+1) = s(end);
end

% figure
% hold on
% plot(time,s)
% plot(time,v)
% plot(time,a)
% plot(time,b)
% xlabel('time')
% legend('s(t)','v(t)','a(t)','b(t)')
% 
% figure
% hold on
% plot(time(1:end-1),diff(s)./diff(time))
% plot(time(1:end-1),diff(v)./diff(time))
% plot(time(1:end-1),diff(a)./diff(time))
% xlabel('time')
% legend('diff: s(t)','diff: v(t)','diff: a(t)')

pos = [time;s];
vel = [time;v];
accel = [time;a];
jerk = [time;b];

