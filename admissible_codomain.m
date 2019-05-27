function y = admissible_codomain(BM,AM,DM,VM,Smin,Smax,s0,v0,a0,sf,vf,af)
%test pripustnosti x0 a xf

%1) otestujem AM a DM
y = 0;
if (a0 <= AM) && (a0 >= -DM) && (af <= AM) && (af >= -DM)
   
else
    'spatne vstupni hodnoty a0 nebo af'
    y = 1;
end;

%2) otestujem VM
if y == 0
    pom1 = a0^2/(2*BM);
    pom2 = af^2/(2*BM);
    vrch1 = v0-pom1;    %VB+
    vrch2 = v0+pom1;    %VB-
    vrch3 = vf-pom2;    %VE+
    vrch4 = vf+pom2;    %VE-
    if (abs(vrch1)<=VM) && (abs(vrch2)<=VM) && (abs(vrch3)<=VM) && (abs(vrch4)<=VM)
        
    else
        'spatne vstupni hodnoty v0 nebo vf'
        y = 1;    
    end;
end;

%otestujem Smin, Smax (krome codomain v-a otestujem zda jsou vstupni polohy
%                        v rozsahu smin, Smax)
if y == 0;
    if (s0 <= Smax) && (s0 >= Smin) && (sf <= Smax) && (sf >= Smin)
   
    else
        'spatne vstupni hodnoty s0 nebo sf'
        y = 1;
    end;
end;