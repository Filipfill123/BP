function y = s_minmax_values(sol,u1,s0,v0,a0,sf,vf,af)
s1_p = s0;
s2_p = sf;
s3_p = s0;                              %zatim dame neutralni hodnotu
m = 0;                                  %rika kolik mame novych s3_p
b = [u1,0,-u1,0,-u1,0,u1,0];
a_0 = a0;
v_0 = v0;
s_0 = s0;
coeff = zeros(8,4);                     %matice koeficientu viz. zprava

% hledame extremy funkce s(t) => musime najit takove casy, kde v(t) = 0
%projizdime vsechny intervaly t1..t7 + jeden krok abych pozdeji mohl zkontrolovat jestli reseni vede do xf.
for k = 1:8                            
    ti = roots([b(k)/2 a_0 v_0]);       %urcim koreny polynomu v(t) v danem casovem intervalu ti
    for i = 1:length(ti)                %pro kazdy REALNY KLADNY koren-
        if (isreal(ti(i))) && (ti(i)>=0) && (ti(i) <= sol(k))%a mensi nez dany casovy interval z reseni 'sol'-
            m = m + 1;
            s3_p(m) = s_0 + v_0*ti(i) + a_0/2*ti(i)^2 + b(k)/6*ti(i)^3; %urcim novy mozny extrem
        end;
    end;
    coeff(k,:) = [b(k) a_0 v_0 s_0];    %ulozim pocatecni hodnoty stavu jako koeff polynomu
    s_0 = s_0 + v_0*sol(k) + a_0/2*sol(k)^2 + b(k)/6*sol(k)^3;%spoctu nove poc hodnoty pro dalsi krok
    v_0 = v_0 + a_0*sol(k) + b(k)/2*sol(k)^2;
    a_0 = a_0 + b(k)*sol(k);
end;
y.smax = max([s1_p,s2_p,s3_p]);
y.smin = min([s1_p,s2_p,s3_p]);
y.coeff = coeff;
%!!!!!!!!!!!!!!! end of s-values
