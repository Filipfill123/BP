function logQ = quatLog(Q)
% Q = [cos(theta/2);sin(theta/2)*dirVect];
% log(Q) = [0;theta/2*dirVect];

theta = 2*acos(Q(1));

if abs(theta) < 1e-10
    logQ = [0;0;0;0];
else
    dirVect = Q(2:4)/sin(theta/2);
    logQ = [0;theta/2*dirVect];
end