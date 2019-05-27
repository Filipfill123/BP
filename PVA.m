function T = PVA(s, vmax, amax)
% generator casovych okamziku prepnuti zrychleni mezi hodnotami -maxAcc, 0,
% +maxAcc
s = abs(s);
if vmax/amax<sqrt(s/amax)
    
t1 = vmax/amax;
t2 = s/vmax;
t3 = vmax/amax + s/vmax;

else
    t1 = sqrt(s/amax);
    t2 = t1;
    t3 = 2*t1;
end

T = [t1;t2;t3;s];