function Quaternion = rotMatrixAndAngularVel2Quaternion(rotMangVel)
% prevod [R,omega,domega] -> [Quat,dQuat,ddQuat]

% polohy
R = rotMangVel(:,1:3);

Quat = [0.5*sqrt(R(1,1) + R(2,2) + R(3,3) + 1);...
        0.5*sign(R(3,2) - R(2,3))*sqrt(abs(R(1,1) - R(2,2) - R(3,3) + 1));...
        0.5*sign(R(1,3) - R(3,1))*sqrt(abs(R(2,2) - R(3,3) - R(1,1) + 1));...
        0.5*sign(R(2,1) - R(1,2))*sqrt(abs(R(3,3) - R(1,1) - R(2,2) + 1))];
Quaternion = Quat;

% rychlosti
if size(rotMangVel,2) > 3
    omega = rotMangVel(:,4);
    
    OMEGA = 0.5*[-Quat(2),-Quat(3),-Quat(4);Quat(1),Quat(4),-Quat(3);-Quat(4),Quat(1),Quat(2);Quat(3),-Quat(2),Quat(1)];
    dQuat = OMEGA*omega;

    Quaternion = [Quaternion,dQuat];
    
    if size(rotMangVel,2) > 4
        domega = rotMangVel(:,5);
        
        dOMEGA = 0.5*[-dQuat(2),-dQuat(3),-dQuat(4);dQuat(1),dQuat(4),-dQuat(3);-dQuat(4),dQuat(1),dQuat(2);dQuat(3),-dQuat(2),dQuat(1)];
        ddQuad = dOMEGA*omega + OMEGA*domega;
        
        Quaternion = [Quaternion,ddQuad];
    end
end

   
