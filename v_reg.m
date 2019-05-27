function [y] = v_reg(Bm,Am,Dm,Vm,x0,xf,t0) %x0 = [s0,v0,a0]  xf = [sf,vf,af]

    vrch(1)=(x0(2)-x0(3)^2/(2*Bm)); %VB+
    vrch(2)=(x0(2)+x0(3)^2/(2*Bm)); %VB-
    vrch(3)=(xf(2)-xf(3)^2/(2*Bm)); %VE+
    vrch(4)=(xf(2)+xf(3)^2/(2*Bm)); %VE-

%3)urcime prvni rizeni (+BM nebo -BM)
%zalezi na sign(a0)
%pro a0>=0: urcime kde lezi End, pokud je pod B-esickem zaciname ridit +, jinak -
if x0(3) >= 0
    if (vrch(4) <= vrch(2))
        u1 = '-BM';
    elseif (vrch(3) < vrch(1)) && (xf(3) > 0)
        u1 = '-BM';
    else
        u1 = 'BM';
    end;
else
%pro a0 < 0: urcime kde lezi End, pokud je pod B-esickem zaciname ridit +, jinak -
    if (vrch(3) >= vrch(1))
        u1 = 'BM';
    elseif (vrch(4) > vrch(2)) && (xf(3) < 0)
        u1 = 'BM';
    else
        u1 = '-BM';
    end;
end;

%protoze x0 a xf prosli kontrolou existence reseni a urcil jsem pocatecni rizeni
%musi v nasledujicich dvou soustavach pro (+BM) nebo (-BM) existovat reseni respektujici omezeni.
%A to jedinne!!???
%spoctu soustavu 1) 2) 3b) (respektuje omezeni), pokud je reseni kladne => OK
%pokud ne, spoctu soustavu 1) 2) 3a) (nerespektuje omezeni), ta ma obecne
%dve reseni (kvadratic rce), z nich podle predpokladu musi byt jedno kladne
%a dokonce respektujici omezeni, protoze prvni soustava kladne reseni nemela
switch u1
    case 'BM'       
        rizeni = [Bm 0 -Bm 0];
% rce 1) 2) 3b)
        sol1 = [(Am - x0(3))/Bm, (2*Bm*(xf(2)-x0(2))+x0(3)^2+xf(3)^2 - 2*Am^2)/(2*Bm*Am), (Am - xf(3))/Bm]; 
        if sol1 >= 0
            intervaly = sol1;
        else
% rce 1) 2) 3a)
            tt3a = (-xf(3)-1/2*(2*xf(3)^2+2*x0(3)^2+4*Bm*xf(2)-4*Bm*x0(2))^(1/2))/Bm;
            tt1a = (xf(3)-x0(3)+Bm*tt3a)/Bm;
            sol2 = [tt1a,0,tt3a];
            if sol2 >= 0
                intervaly = sol2;
            else
                tt3b = (-xf(3)+1/2*(2*xf(3)^2+2*x0(3)^2+4*Bm*xf(2)-4*Bm*x0(2))^(1/2))/Bm;
                tt1b = (xf(3)-x0(3)+Bm*tt3b)/Bm;
                intervaly = [tt1b,0,tt3b];
            end;
        end;
t1 =intervaly(1);
t2 =intervaly(2);
t3 =intervaly(3);
sr = x0(2)*(t1+t2+t3)+(t1*t3+t2*t3+1/2*t2^2+1/2*t1^2+t1*t2+1/2*t3^2)*x0(3)+1/2*Bm*t1*t3^2+1/2*Bm*t1^2*t3+x0(1)-1/6*Bm*t3^3+1/2*Bm*t1*t2^2+1/6*Bm*t1^3+1/2*Bm*t1^2*t2+Bm*t1*t2*t3;        

%to same pro -BM. Stejny pristup s podobnymi rovnicemi pro rizeni [-Bm 0 Bm 0]       
    case '-BM'
        rizeni = [-Bm 0 Bm 0];
% rce 1) 2) 3b)
        sol1 = [(x0(3)-Dm)/Bm, (2*Bm*(xf(2)-x0(2))-x0(3)^2-xf(3)^2 + 2*Dm^2)/(2*Bm*Dm), (xf(3)-Dm)/Bm];
        if sol1 >= 0
            intervaly = sol1;
        else        
% rce 1) 2) 3a)
            tt3a = (xf(3)-1/2*(2*xf(3)^2+2*x0(3)^2+4*Bm*x0(2)-4*Bm*xf(2))^(1/2))/Bm;
            tt1a = (x0(3)-xf(3)+Bm*tt3a)/Bm;
            sol2 = [tt1a,0,tt3a];
            if sol2 >= 0
                intervaly = sol2;
            else
            tt3b = (xf(3)+1/2*(2*xf(3)^2+2*x0(3)^2+4*Bm*x0(2)-4*Bm*xf(2))^(1/2))/Bm;
            tt1b = (x0(3)-xf(3)+Bm*tt3b)/Bm;
            intervaly = [tt1b,0,tt3b];
            end;
        end;      
t1 =intervaly(1);
t2 =intervaly(2);
t3 =intervaly(3);
sr = x0(2)*(t1+t2+t3)+(1/2*t1^2+t1*t3+t2*t3+t1*t2+1/2*t3^2+1/2*t2^2)*x0(3)-1/6*Bm*t1^3-Bm*t1*t2*t3-1/2*Bm*t1^2*t3+x0(1)+1/6*Bm*t3^3-1/2*Bm*t1*t3^2-1/2*Bm*t1*t2^2-1/2*Bm*t1^2*t2;

end;




casy(1) = t0;
casy(2) = intervaly(1) + casy(1);
casy(3) = intervaly(2) + casy(2);
casy(4) = intervaly(3) + casy(3);
y.int = intervaly;
y.riz = rizeni;
y.cas = casy;
y.sr = sr;