function T = BAVS(S,v_max,a_max,b_max)

jerk = b_max;
c_acc = a_max;
c_vel = v_max;
c_dec = a_max;
c_pos_p = inf;
c_pos_m = -inf;
s_init = 0;
v_init = 0;
a_init = 0;
s_fin = S;
v_fin = 0;
a_fin = 0;%sqrt(a_init^2+2*jerk*(v_fin-v_init))+0;
% c = 1
% s = s_init + v_init*c + a_init/2*c^2 + jerk/6*c^3
% v = v_init + a_init*c + jerk/2*c^2
% a = a_init + jerk*c
% return

%t-opt regulator rychlosti => neni potreba pro samotny algoritmus
av_reg = v_reg(jerk,c_acc,c_dec,c_vel,[s_init,v_init,a_init],[s_fin,v_fin,a_fin],0); %x0 = [s0,v0,a0]  xf = [sf,vf,af]
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
%return

tic
%A)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%1) test of inputs = zda je x0 a xf v pripustne oblasti
error1 = admissible_codomain(jerk,c_acc,c_dec,c_vel,c_pos_m,c_pos_p,s_init,v_init,a_init,s_fin,v_fin,a_fin);
if error1 == 1
    'error1'
    return
end;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of 1)

%2) solve U1+
U1_plus = U1(jerk,c_acc,c_dec,c_vel,s_init,v_init,a_init,s_fin,v_fin,a_fin);
if U1_plus.numsol > 0                   %pokud existuje real kladne res
    for i = 1:U1_plus.numsol            %projdi je 
        tau7(i) = sum(U1_plus.sol(i,:));%a urci tau7 = t_final
    end;
    [t_minU1p,indexU1p] = min(tau7);    %vyber minimalni (uloz cas a index reseni)
    solU1p = U1_plus.sol(indexU1p,:);   %uloz minimalni reseni vypoctene z U1+
    U1p_mares = true;                   %nejake reseni zde bylo nalezeno
else
    U1p_mares = false;                  %reseni v U1+ neni
end;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of 2)

%3) change of input variables + solve U1-
%pozor jsou zde zamenena znaminka a poradi nekterych vstupnich parametru, porovnej s U1+
U1_minus = U1(jerk,c_dec,c_acc,c_vel,-s_init,-v_init,-a_init,-s_fin,-v_fin,-a_fin);
if U1_minus.numsol > 0                  %stejne jako u U1+
    for i = 1:U1_minus.numsol
        tau7m = sum(U1_minus.sol(i,:));
    end;
    [t_minU1m,indexU1m] = min(tau7m);
    solU1m = U1_minus.sol(indexU1m,:);
    U1m_mares = true;
else
    U1m_mares = false;
end;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of 3)

%4)!!!!!!!!!!!!!!!!!!!choose t-opt sol
if U1p_mares && U1m_mares                   %kdyz je v U1+ i U1- reseni
    if t_minU1m < t_minU1p                  %vyberu rychlejsi
        solution = solU1m;
        u = [-jerk,0,jerk,0,jerk,0,-jerk,0];%priradim prislusnou posloupnost rizeni
    else
        solution = solU1p;
        u = [jerk,0,-jerk,0,-jerk,0,jerk,0];  
    end;
elseif U1p_mares                            %kdyz je reseni jen v U1+
    solution = solU1p;                      %vyberu ho
    u = [jerk,0,-jerk,0,-jerk,0,jerk,0];  
elseif U1m_mares                            %naopak pro U1-
    solution = solU1m;
    u = [-jerk,0,jerk,0,jerk,0,-jerk,0];
else                                        %reseni neni nikde, CHYBA ALGORITMU, protoze z teorie -
    'error6 - unexpected'                   %reseni musi existovat pro pripustne vstupy a ty uz byly - 
    return                                  %zkontrolovany. Nastalo pri chybe malych cisel => osetreno, ale...:-(
end;
%!!!!!!!!!!!!!!!!!!! end of t_opt

%5)!!!!!!!!!!!!!!!!! min-max values
a_minmax = a_minmax_values(solution,u(1),a_init,a_fin);                              %urci a_min,a_max
v_minmax = v_minmax_values(solution,u(1),v_init,a_init,v_fin,a_fin);                 %urci v_min,v_max

%urci s_min,s_max,a matici koeficientu polynomu (vice ve zprave) -
%POZOR prvni parametr solution je prodlouzen
s_minmax = s_minmax_values([solution 0],u(1),s_init,v_init,a_init,s_fin,v_fin,a_fin);
%!!!!!!!!!!!!! end of min-max values

%6)kontrola zda reseni opravdu vede do cile (do xf)
tol = 10^-7;                                   %tolerance pro numer chyby
kontrola = s_minmax.coeff(8,:);
if (abs(kontrola(2) - a_fin) < tol) && (abs(kontrola(3) - v_fin) < tol) && (abs(kontrola(4) - s_fin) < tol)
    
else
    'error 7 - chyba malych cisel'
    return
end;
%6)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%7)!!!!!!znovu kontrola omezeni VM,AM,DM pro vysledne reseni
%melo by to vzdy platit, ale...
if (a_minmax.amax <= c_acc+tol) && (a_minmax.amin >= -c_dec-tol) && (v_minmax.vmax <= c_vel+tol) && (v_minmax.vmin >= -c_vel-tol)

else                                        %kontrola omezeni AM,VM,DM neprosla
    'error 6 - unexpected - jak je to mozny'%CHYBA ALGORITMU
    return
end;

%8)!!!!!!!!!!!!!! test c_pos
%kontrola zda se reseni "vejde" do pracovniho prostoru S
if (s_minmax.smin >= c_pos_m) && (s_minmax.smax <= c_pos_p)
    for i = 1:7                             %vypocet vektoru tau (tau1...tau7)
        tau(i) = sum(solution(1:i));
    end;
 
%     'NALEZENE RESENI'
%     solution                                %vypis reseni
%     tau                                     %vypis tau
%     u                                       %vypis prislusneho rizeni
%     finaltime = tau(7)                      %vypis vysledneho casu
%     a_minmax                                %vypis min-max hodnot
%     v_minmax
%     s_minmax
        
else        %nalezene reseni vede do cile ale nevejde se do S
    'error 5 - poruseni c_pos'
    s_minmax
    return
end;

T = tau;
