function y = a_minmax_values(sol,u1,a0,af)
%urceni minimalnich a maximalnich hodnot t-opt reseni;
a1_p = a0;
a2_p = af;
a3_p = a0 + u1*sol(1);
a4_p = af - u1*sol(7);
y.amin = min([a1_p,a2_p,a3_p,a4_p]);
y.amax = max([a1_p,a2_p,a3_p,a4_p]);
%!!!!!!!!!!!!!!!!!! end of a-values
