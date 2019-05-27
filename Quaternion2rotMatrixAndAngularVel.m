function rotMangVel = Quaternion2rotMatrixAndAngularVel(Quaternion)
% prevod [Quat,dQuat,ddQuat] -> [R,omega,domega]

% polohy
Quat = Quaternion(:,1);

R = [2*(Quat(1)^2 + Quat(2)^2) - 1,         2*(Quat(2)*Quat(3) - Quat(1)*Quat(4)), 2*(Quat(2)*Quat(4) + Quat(1)*Quat(3));...
     2*(Quat(2)*Quat(3) + Quat(1)*Quat(4)), 2*(Quat(1)^2 + Quat(3)^2) - 1,         2*(Quat(3)*Quat(4) - Quat(1)*Quat(2));...
     2*(Quat(2)*Quat(4) - Quat(1)*Quat(3)), 2*(Quat(3)*Quat(4) + Quat(1)*Quat(2)), 2*(Quat(1)^2 + Quat(4)^2) - 1];
rotMangVel = R;

% rychlosti
if size(Quaternion,2) > 1
    dQuat = Quaternion(:,2);
    
   OMEGA = 0.5*[-Quat(2),-Quat(3),-Quat(4);Quat(1),Quat(4),-Quat(3);-Quat(4),Quat(1),Quat(2);Quat(3),-Quat(2),Quat(1)];
    omega = 4*OMEGA'*dQuat;
    
    rotMangVel = [rotMangVel,omega];
    
    if size(Quaternion,2) > 2
        ddQuat = Quaternion(:,3);
        
        dOMEGA = 0.5*[-dQuat(2),-dQuat(3),-dQuat(4);dQuat(1),dQuat(4),-dQuat(3);-dQuat(4),dQuat(1),dQuat(2);dQuat(3),-dQuat(2),dQuat(1)];
        domega = 4*OMEGA'*(ddQuat - dOMEGA*omega);
        
        rotMangVel = [rotMangVel,domega];
    end
end
   
