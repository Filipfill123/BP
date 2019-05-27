clear
clc
close all

S = 1.1;
maxVelocity = 0.5;
maxAcceleration = 1;
res = 0.17;

[pos,vel,accel] = PVA_genTrajectories_rounding(S,maxVelocity,maxAcceleration,res)

figure
hold on
plot(pos(1,:),pos(2,:),'b')
plot(pos(1,:),vel(2,:),'g')
plot(pos(1,:),accel(2,:),'r')