function [pos,vel,accel] = PVA_genTrajectories(S,maxVelocity,maxAcceleration,res)

% [pos,vel,accel] = PVA_genTrajectories(S,maxVelocity,maxAcceleration,res)
%
% generator primkoveho pohybu po draze (S) s omezenim na rychlost
% (maxVelocity) a zrychleni (maxAcceleration)
% res =  interval mezi casovymmi vzorky (priblizne)

maxVelocity = abs(maxVelocity);
maxAcceleration = abs(maxAcceleration);

T = PVA(S,maxVelocity,maxAcceleration);

if T(1)~=T(2)
    int1 = linspace(0,T(1),round(T(1)/res));
    int2 = linspace(T(1),T(2),round((T(2)-T(1))/res));
    int3 = linspace(T(2),T(3),round((T(3)-T(2))/res));
    int2 = int2(2:end-1);
else
    int1 = linspace(0,T(1),round(T(1)/res));
    int2 = [];
    int3 = linspace(T(2),T(3),round((T(3)-T(2))/res));
    int3 = int3(2:end);
end



accel = [int1,int2,int3;...
         maxAcceleration*ones(1,length(int1)),0*ones(1,length(int2)),-maxAcceleration*ones(1,length(int3))];
vel = [int1,int2,int3;...
       maxAcceleration*int1,maxAcceleration*T(1)*ones(1,length(int2)),-maxAcceleration*(int3-T(2)-T(1))];
pos = [int1,int2,int3;...
         0.5*maxAcceleration*int1.^2,maxAcceleration*(T(1)*int2-0.5*T(1)^2),-0.5*maxAcceleration*(int3.^2-T(2)^2)+maxAcceleration*(T(1)+T(2))*(int3-T(2))+maxAcceleration*(T(1)*T(2)-0.5*T(1)^2)];
     
pos(2,:) = sign(S)*pos(2,:);
vel(2,:) = sign(S)*vel(2,:);
accel(2,:) = sign(S)*accel(2,:);