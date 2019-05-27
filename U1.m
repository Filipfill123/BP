function y = U1(BM,AM,DM,VM,s0,v0,a0,sf,vf,af)
e = 10^-13;                 %tolerance pro chyby malych cisel
A = zeros(1,7);             
Apositiv = zeros(1,7);
solution = zeros(1,7);
index = 0;
pocet = 0;
num_of_sol = 0;

%set1
s1t1 = (AM - a0)/BM;
s1t2 = (2*BM*(VM-v0)+a0^2-2*AM^2)/2/BM/AM;
s1t3 = AM/BM;
s1t4 = 1/24*(12*((-v0-VM)*BM+1/2*a0^2)*DM*AM^2+(((-vf-VM)*12*BM+6*af^2)*DM^2+((sf-s0)*24*BM^2+(v0*a0-vf*af)*24*BM+8*(af^3-a0^3))*DM+(vf^2-VM^2)*12*BM^2-12*af^2*BM*vf+3*af^4)*AM+((v0^2-VM^2)*12*BM^2-12*a0^2*BM*v0+3*a0^4)*DM)/DM/AM/BM^2/VM;   
s1t5 = DM/BM;
s1t6 = (2*BM*(VM-vf)+af^2-2*DM^2)/2/BM/DM;
s1t7 = (af+DM)/BM;
index = index+1;
A(index,:) = [s1t1 s1t2 s1t3 s1t4 s1t5 s1t6 s1t7];

if A(index,:) >= 0
    pocet = pocet+1;
    Apositiv(pocet,:) = A(index,:);
    num_of_sol = num_of_sol+1;
    solution(num_of_sol,:) = A(index,:);
end;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%set2
s2t1 = s1t1;
s2t2 = s1t2;
s2t3 = s1t3;
s2t4 = 1/24*(6*(2*(-vf-VM)*BM+af^2)*AM*(2*af^2+4*BM*(VM-vf))^(1/2)+6*(2*(-VM-v0)*BM+a0^2)*AM^2+8*(3*(sf-s0)*BM^2+3*(-af*vf+a0*v0)*BM+af^3-a0^3)*AM+12*(v0^2-VM^2)*BM^2-12*BM*v0*a0^2+3*a0^4)/AM/BM^2/VM;
s2t5 = (2*af^2+4*BM*(VM-vf))^(1/2)/BM/2;
s2t6 = 0;
s2t7 = af/BM + s2t5;

index = index+1;
A(index,:) = [s2t1 s2t2 s2t3 s2t4 s2t5 s2t6 s2t7];

if A(index,:) >= 0
    pocet  = pocet+1;
    Apositiv(pocet,:) = A(index,:);
    ok = testDM(BM,DM,af,s2t7,e);
    if ok == 1
        num_of_sol = num_of_sol+1;
        solution(num_of_sol,:) = A(index,:);
    else
%        'set2_1 je kladny ale necti DM'
    end;   
end;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%set3
s3t1 = s1t1;
s3t3 = s1t3;
s3t4 = 0;
s3t5 = s1t5;
s3t7 = s1t7;

p0 = 12*(2*(s0-sf)*AM-v0^2+vf^2)*BM^2+12*((v0+vf)*AM^2+2*(-a0*v0+vf*af+2*vf*DM)*AM+a0^2*v0-af^2*vf+2*DM^2*vf)*BM+6*(-a0^2+2*DM^2-af^2)*AM^2+8*(-af^3+3*DM^3+a0^3-3*af^2*DM)*AM+3*af^4+12*DM^4-12*DM^2*af^2-3*a0^4;
p1 = 24*(DM^2+1/2*AM*DM+vf*BM-1/2*af^2)*(AM+DM)*BM;
p2 = 12*DM*BM^2*(AM+DM);
p = [p2 p1 p0];
% for i = 1:3
%     if abs(p(i)) < e
%         p(i) = 0;
%     end;
% end;
p = facelift(p);
r3y = roots(p);

for i = 1:length(r3y)
    s3t6 = r3y(i);
    s3t2 = (2*BM*(s3t6*DM-v0+vf)+a0^2-af^2+2*(DM^2-AM^2))/2/BM/AM;
    A(index+i,:) = [s3t1 s3t2 s3t3 s3t4 s3t5 s3t6 s3t7];
    
    if A(index+i,:) >= 0 & isreal(A(index+i,:))
        pocet = pocet+1;
        Apositiv(pocet,:) = A(index+i,:);
        ok = v_minmax_values(A(index+i,:),BM,v0,a0,vf,af);
        if ok.vmax <= VM+e && ok.vmin >= -VM-e
            num_of_sol = num_of_sol+1;
            solution(num_of_sol,:) = A(index+i,:);
        else
%            'set3_i je kladny ale necti VM'
        end;   
    end;
end;
index = index+length(r3y);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%set4
s4t1 = s1t1;
s4t4 = 0;
s4t5 = 0;
s4t6 = 0;

p0 = ((v0+vf)*2*BM+af^2-a0^2)*6*AM^2+((s0-sf)*3*BM^2+(-af*vf-a0*v0)*3*BM-af^3+a0^3)*8*AM+(-v0^2+vf^2)*12*BM^2+(af^2*vf+a0^2*v0)*12*BM+3*af^4-3*a0^4;
p1 = -48*BM*(1/2*af^2-1/2*AM*af+BM*vf)*(af-AM);
p2 = 12*BM^2*(-6*AM*af+AM^2+5*af^2+2*BM*vf);
p3 = -24*BM^3*(-AM+2*af);
p4 = 12*BM^4;
p = [p4 p3 p2 p1 p0];
% for i = 1:5
%     if abs(p(i)) < e
%         p(i) = 0;
%     end;
% end;
p = facelift(p);
r4y = roots(p);

for i = 1:length(r4y)
    s4t7 = r4y(i);
    s4t3 = s1t3 + s4t7 - af/BM;
    s4t2 = (2*BM^2*s4t7^2+2*BM*(vf-v0-2*s4t7*af)+a0^2+af^2-2*AM^2)/2/BM/AM;
    A(index+i,:) = [s4t1 s4t2 s4t3 s4t4 s4t5 s4t6 s4t7]; 
    
    if A(index+i,:) >= 0 & isreal(A(index+i,:))
        pocet = pocet+1;
        Apositiv(pocet,:) = A(index+i,:);
        ok1 = v_minmax_values(A(index+i,:),BM,v0,a0,vf,af);
        ok2 = testDM(BM,DM,af,s4t7,e);
        if ok1.vmax <= VM+e && ok1.vmin >= -VM-e && ok2 == 1
            num_of_sol = num_of_sol+1;
            solution(num_of_sol,:) = A(index+i,:);
        else
%            'set4_i je kladny ale necti VM nebo DM'
        end;   
    end;
end;
index = index+length(r4y);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%set5
s5t3 = ((2*a0^2+4*BM*(VM-v0))^(1/2))/2/BM;
s5t1 = -a0/BM + s5t3;
s5t2 = 0;
s5t4 = (6*(2*(-v0-VM)*BM+a0^2)*DM*(2*a0^2+4*BM*(VM-v0))^(1/2)+6*(2*(-vf-VM)*BM+af^2)*DM^2+8*(3*(-s0+sf)*BM^2+3*(v0*a0-vf*af)*BM+af^3-a0^3)*DM+12*(vf^2-VM^2)*BM^2-12*BM*vf*af^2+3*af^4)/DM/BM^2/VM/24;
s5t5 = s1t5;
s5t6 = s1t6;
s5t7 = s1t7;
index = index+1;
A(index,:) = [s5t1 s5t2 s5t3 s5t4 s5t5 s5t6 s5t7];

if A(index,:) >= 0
    pocet = pocet+1;
    Apositiv(pocet,:) = A(index,:);
    ok = testAM(BM,AM,a0,s5t1,e);
    if ok == 1
        num_of_sol = num_of_sol+1;
        solution(num_of_sol,:) = A(index,:);
    else
%        'set5 je sice kladny ale necti AM'
    end;
end;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%set6
t1 = s5t1;
t2 = 0;
t3 = s5t3;
t5 = s2t5;
t6 = 0;
t7 = s2t7;
t4 = (2*(3*af^2-6*BM*(vf+VM))*t5*BM+2*(3*a0^2-6*BM*(v0+VM))*t3*BM+4*af^3-12*BM*vf*af-4*a0^3+12*a0*BM*v0-12*BM^2*(s0-sf))/BM^2/VM/12;
A(index+i,:) = [t1 t2 t3 t4 t5 t6 t7]; 
    
if A(index+i,:) >= 0
    pocet = pocet+1;
    Apositiv(pocet,:) = A(index+i,:);
    ok1 = v_minmax_values(A(index+i,:),BM,v0,a0,vf,af);
    ok2 = testDM(BM,DM,af,t7,e);
    ok3 = testAM(BM,AM,a0,t1,e);
    if ok1.vmax <= VM+e && ok1.vmin >= -VM-e && ok2 == 1 && ok3 == 1
        num_of_sol = num_of_sol+1;
        solution(num_of_sol,:) = A(index+i,:);
    else
%       'set6_i je kladny ale necti VM nebo DM nebo AM'
    end;   
end;
index = index+4;
%!!!!!!!!!!!!!!!!!!!!!!!!

%set7 zkontrolovat AM,VM
s7t2 = 0;
s7t4 = 0;
s7t5 = 0;
s7t7 = s1t7;

p0 = ((v0+vf)*2*BM+a0^2-af^2)*6*DM^2+((s0-sf)*3*BM^2+(a0*v0+vf*af)*3*BM+a0^3-af^3)*8*DM+(v0^2-vf^2)*12*BM^2+(vf*af^2+v0*a0^2)*12*BM-3*af^4+3*a0^4;
p1 = 48*(DM+a0)*(1/2*DM*a0+BM*v0+1/2*a0^2)*BM;
p2 = 12*BM^2*(DM^2+6*DM*a0+5*a0^2+2*BM*v0);
p3 = 24*BM^3*(2*a0+DM);
p4 = 12*BM^4;
p = [p4 p3 p2 p1 p0];
% for i = 1:5
%     if abs(p(i)) < e
%         p(i) = 0;
%     end;
% end;
p = facelift(p);
r7y = roots(p);

for i = 1:length(r7y)
    s7t1 = r7y(i);
    s7t3 = DM/BM + a0/BM + s7t1;
    s7t6 = (2*BM^2*s7t1^2+2*BM*(v0-vf+2*a0*s7t1)+a0^2+af^2-2*DM^2)/2/BM/DM; 
    A(index+i,:) = [s7t1 s7t2 s7t3 s7t4 s7t5 s7t6 s7t7];

    if A(index+i,:) >= 0 & isreal(A(index+i,:))
        pocet = pocet+1;
        Apositiv(pocet,:) = A(index+i,:);
        ok1 = v_minmax_values(A(index+i,:),BM,v0,a0,vf,af);
        ok2 = testAM(BM,AM,a0,s7t1,e);
        if ok1.vmax <= VM+e && ok1.vmin >= -VM-e && ok2 == 1
            num_of_sol = num_of_sol+1;
            solution(num_of_sol,:) = A(index+i,:);
        else
%            'set7_i je kladny ale necti VM nebo AM'
        end;   
    end;
end;
index = index+length(r7y);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%set8
s8t2 = 0;
s8t4 = 0;
s8t5 = 0;
s8t6 = 0;

p0 = 36*af^2*a0^2*BM*v0+48*af^3*BM*v0*a0-144*BM^3*v0*a0*s0-72*BM^2*vf*a0^2*v0-36*a0^2*af^2*BM*vf-144*BM^3*vf*af*sf+144*BM^3*vf*af*s0+48*BM*vf*af*a0^3-a0^6+17*af^6+9*a0^2*af^4+48*af^3*BM^2*sf-18*BM*v0*af^4+72*BM^2*v0*af^2*vf+72*BM^4*s0^2+72*BM^3*v0^3-72*BM^3*vf^3+72*BM^4*sf^2-9*af^2*a0^4+6*BM*v0*a0^4-144*BM^4*s0*sf-72*BM^3*v0*vf^2+72*BM^3*vf*v0^2+180*BM^2*vf^2*af^2-48*a0^3*BM^2*sf-48*af^3*BM^2*s0-102*BM*vf*af^4-36*BM^2*v0^2*a0^2+48*a0^3*BM^2*s0+36*a0^2*BM^2*vf^2-36*af^2*BM^2*v0^2-144*BM^2*vf*af*v0*a0-16*af^3*a0^3+144*BM^3*v0*a0*sf+18*BM*vf*a0^4;
p1 = (288*BM^3*vf*s0+288*BM^2*vf^2*af-48*af^2*a0^3-288*BM^3*vf*sf-288*BM^2*vf*v0*a0+144*af^2*BM*v0*a0-240*af^3*BM*vf+48*af^5+144*af^2*BM^2*sf-144*af^2*BM^2*s0+96*BM*vf*a0^3);
p2 = (72*BM^2*v0^2+72*BM^2*vf^2-36*af^2*a0^2+18*a0^4+72*a0^2*BM*vf+72*af^2*BM*v0-72*a0^2*BM*v0-144*BM^2*v0*vf-72*af^2*BM*vf+18*af^4);
p3 = (-144*BM^2*sf-144*BM*v0*a0-48*af^3+48*a0^3+144*BM^2*s0+144*BM*vf*af);
p4 = (72*BM*vf-72*BM*v0-36*af^2+36*a0^2);
p = [p4 p3 p2 p1 p0];
% for i = 1:5
%     if abs(p(i)) < e
%         p(i) = 0;
%     end;
% end;
p = facelift(p);
r8dm = roots(p);

for i = 1:length(r8dm)
    dm = r8dm(i);
    am = -((144*s0*sf*dm-72*s0^2*dm-72*sf^2*dm+144*vf^2*s0-144*vf*v0*s0+144*vf*v0*sf-144*vf^2*sf)*BM^4+(144*vf*af*sf*dm-72*vf*dm^2*sf-144*sf*v0*a0*dm+72*vf^3*dm+72*v0*af^2*s0-72*v0*dm^2*s0+72*vf*dm^2*s0+144*vf^3*af-144*vf^2*v0*af-72*vf*a0^2*sf-216*vf^2*dm*v0+144*s0*v0*a0*dm+216*vf*dm*v0^2+72*vf*a0^2*s0-72*v0*af^2*sf+144*vf*af^2*sf-72*v0^3*dm-144*vf^2*v0*a0+144*vf*v0^2*a0+72*v0*dm^2*sf-144*vf*af^2*s0-144*vf*af*s0*dm)*BM^3+(-48*a0^3*s0*dm+48*a0^3*sf*dm-36*a0^2*af^2*s0+36*v0^2*dm*a0^2-108*v0^2*dm*af^2-36*a0^2*dm^2*sf+36*a0^2*dm^2*s0+36*sf*a0^2*af^2-180*vf^2*dm*af^2+108*vf^2*dm*a0^2-144*vf*dm^3*v0-72*v0^2*a0*af^2-120*v0*a0^3*vf+72*vf^2*dm^2*af-216*vf*dm*v0*a0^2+216*vf*dm*v0*af^2-72*v0*dm^2*vf*af-72*vf*dm^2*v0*a0+72*v0^2*dm^2*a0+144*v0*a0*vf*af^2+144*vf*af*v0*a0*dm-192*vf^2*af^3-36*sf*af^4+36*af^4*s0+72*vf^2*dm^3+72*v0^2*dm^3+48*a0^3*vf^2+72*vf^2*af*a0^2+120*vf*af^3*v0+36*af^2*dm^2*sf-36*af^2*dm^2*s0+48*af^3*s0*dm-48*af^3*sf*dm)*BM^2+(72*vf*dm^3*a0^2-72*vf*dm^3*af^2-36*af^4*v0*a0-6*v0*dm*a0^4-54*v0*dm*af^4+24*v0*dm^2*af^3-60*v0*dm^2*a0^3-72*v0*dm^3*a0^2-108*vf*dm*a0^2*af^2+108*v0*dm*a0^2*af^2+36*a0^2*dm^2*vf*af+36*af^2*dm^2*v0*a0-48*vf*af*a0^3*dm-48*af^3*v0*a0*dm+24*vf*a0^5+72*v0*dm^3*af^2+54*vf*dm*a0^4+102*vf*dm*af^4-60*vf*a0^2*af^3-48*vf*af^2*a0^3+60*v0*af^2*a0^3+24*vf*dm^2*a0^3-60*vf*dm^2*af^3+84*vf*af^5-24*v0*af^5)*BM+18*a0^4*dm^3+12*a0^5*dm^2-36*a0^2*dm^3*af^2+12*a0^2*af^5+12*af^5*dm^2-17*af^6*dm+a0^6*dm+18*af^4*dm^3-12*a0^2*dm^2*af^3+27*a0^2*dm*af^4+12*af^4*a0^3-12*af^2*dm^2*a0^3+16*af^3*a0^3*dm-27*a0^4*dm*af^2-12*a0^5*af^2-12*af^7)/((-144*s0*sf+72*s0^2+72*sf^2)*BM^4+(144*vf*af*s0-144*s0*v0*a0+72*v0^3-72*v0*vf^2-144*vf*af*sf+144*sf*v0*a0-72*v0^2*vf+72*vf^3)*BM^3+(-48*a0^3*sf-36*vf^2*af^2-36*v0^2*a0^2+48*a0^3*s0+36*vf^2*a0^2-144*vf*af*v0*a0+36*af^2*v0^2+72*vf*v0*a0^2-48*af^3*s0+48*af^3*sf+72*vf*v0*af^2)*BM^2+(48*vf*af*a0^3+6*vf*af^4-36*af^2*v0*a0^2+6*v0*a0^4+48*af^3*v0*a0-18*af^4*v0-18*vf*a0^4-36*vf*a0^2*af^2)*BM+9*a0^2*af^4+9*a0^4*af^2-16*af^3*a0^3-af^6-a0^6);
    s8t1 = (am - a0)/BM;
    s8t3 = (am + dm)/BM;
    s8t7 = (dm + af)/BM;
    A(index + i,:) = [s8t1 s8t2 s8t3 s8t4 s8t5 s8t6 s8t7];
    if A(index + i,:) >= 0 & isreal(A(index+i,:))
        pocet = pocet+1;
        Apositiv(pocet,:) = A(index+i,:);       
        ok = v_minmax_values(A(index+i,:),BM,v0,a0,vf,af);
        if (ok.vmax <= VM+e) && (ok.vmin >= -VM-e) && (dm <= DM+e) && (am <= AM+e)
            num_of_sol = num_of_sol+1;
            solution(num_of_sol,:) = A(index+i,:);
        else
%           'set8_i je kladny ale necti VM nebo AM nebo DM'
        end;   
    end;
end;
index = index+length(r8dm);

%var8_2
s8t2 = 0;
s8t4 = 0;
s8t5 = 0;
s8t6 = 0;

p0 = -12*a0^2*BM*vf+12*a0^2*BM*v0-3*af^4+12*af^2*BM*vf-12*af^2*BM*v0-12*BM^2*vf^2+24*BM^2*vf*v0-12*BM^2*v0^2+6*a0^2*af^2-3*a0^4;
p1 = -48*BM^3*sf-48*BM^2*v0*a0+48*BM^3*s0-16*BM*af^3+48*BM^2*af*vf+16*BM*a0^3;
p2 = 48*BM^3*v0-24*BM^2*a0^2-24*BM^2*af^2+48*BM^3*vf;
p3 = 0;
p4 = 12*BM^4;
p = [p4 p3 p2 p1 p0];

p = facelift(p);
r8t3 = roots(p);

for i = 1:length(r8t3)
    s8t3 = r8t3(i);
    s8t1 = 1/6*(24*af*BM*vf-12*BM*vf*a0-8*af^3-12*BM*v0*a0+6*BM^3*s8t3^3-9*BM*a0^2*s8t3+30*BM^2*s8t3*vf-15*BM*s8t3*af^2+18*BM^2*v0*s8t3+2*a0^3+24*BM^2*s0-24*BM^2*sf+6*af^2*a0)/BM/(a0^2-af^2+2*BM*(vf-v0));
    s8t7 = -s8t1+s8t3+(af-a0)/BM;
    A(index + i,:) = [s8t1 s8t2 s8t3 s8t4 s8t5 s8t6 s8t7];

    if A(index + i,:) >= 0 & isreal(A(index+i,:))
        pocet = pocet+1;
        Apositiv(pocet,:) = A(index+i,:);       
        ok = v_minmax_values(A(index+i,:),BM,v0,a0,vf,af);
        oka = testAM(BM,AM,a0,s8t1,e);
        okd = testDM(BM,DM,af,s8t7,e);
        if (ok.vmax <= VM+e) && (ok.vmin >= -VM-e) && (okd == 1) && (oka == 1)
            num_of_sol = num_of_sol+1;
            solution(num_of_sol,:) = A(index+i,:);
        else
%            'set8_i je kladny ale necti VM nebo AM nebo DM'
        end;   
    end;
end;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
y.A = A;
y.pos = Apositiv;
y.sol = solution;
y.numsol = num_of_sol;
test = roots([3*BM^3,0, 24*BM^2*v0-12*BM*a0^2, -12*BM*v0*a0-6*a0^2*af+2*af^3+12*BM*v0*af+4*a0^3+12*BM^2*s0-12*BM^2*sf]);
y.test =test;