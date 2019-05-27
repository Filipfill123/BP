%test AM
function y = testAM(bm,am,a0,t1,e)
am_prouzek = a0 + bm*t1;
if am_prouzek <= am + e;
    y = 1;
else
    y = 0;
end;
    