function y = v_minmax_values(sol,u1,v0,a0,vf,af)
%urceni minimalnich a maximalnich hodnot t-opt reseni;

a_2 = a0 + u1*sol(1);   %=a2_p      urceni a(t) v case t2
a_3 = a_2 - u1*sol(3);  %=a3_p      urceni a(t) v case t3
v_2 = v0 + a0*(sol(1)+sol(2)) + u1*(sol(1)^2/2 + sol(1)*sol(2));%urceni v(t) v case t2

%!!!!!!!!!!!!!!!!!! v-values
%!!!!!!!!!v2_p
if sign(a0)~=sign(a_2)
    v2_p = v0 - a0^2/(2*u1);    
else
    v2_p = v0;
end;
%!!!!!!!!!!!!!!!!!!!!

%!!!!!!!!!!!! v3_p
if sign(a_2) ~= sign(a_3)        
    v3_p = v_2 + a_2^2/(2*u1);     
else                               
    v3_p = v0;              
end;                        
%!!!!!!!!!!!!!!!!!!!!!!

%!!!!!!!!!!!! v4_p:   obdobne jako u v0
if sign(af)~=sign(a_3)
    v4_p = vf - af^2/(2*u1);
else
    v4_p = vf;
end;
%!!!!!!!!!!!!!!!!!!

y.vmax = max([v0,v2_p,v3_p,v4_p,vf]);
y.vmin = min([v0,v2_p,v3_p,v4_p,vf]);
%!!!!!!!!!!!!!!!!!! end of v_values