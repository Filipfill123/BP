%test DM
function y = testDM(bm,dm,af,t7,e)
dm_prouzek = af - bm*t7;
if dm_prouzek >= -dm - e;
    y = 1;
else
    y = 0;
end;
    