function expQ = quatExp(Q)
% Q = [0;theta/2*dirVect];
% expQ = [cos(theta/2);sin(theta/2)*dirVect];

theta = 2*norm(Q(2:4));

if abs(theta) < 1e-10
    expQ = [1;0;0;0];
else
    dirVect = Q(2:4)/(theta/2);
    expQ = [cos(theta/2);sin(theta/2)*dirVect];
end