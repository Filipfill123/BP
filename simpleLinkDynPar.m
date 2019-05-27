function [m,T,I] = simpleLinkDynPar(L,d,rho,M)
% rameno je dano jako valcova plna tyc o prumeru d, hustote rho a delce L
% s.s. F ramene je uvazovana s pocatkem na jeho konci s osou x totoznou s
% osou valce!
% 
% L ... delka ramene
% d ... prumer ramene
% rho ... hustota materialu
% M ... hmotnost (motoru - hmotny bod) na konci ramene (v pocatku s.s. ramene)
% 
% m ... hmotnost ramene
% T ... souradnice teziste v s.s. F
% I ... tensor setrvacnosti v tezisti v s.s. F

mL = pi/4*d^2*rho*L;
ILyz = 1/12*mL*(3/4*d^2 + L^2);
ILx = 1/2*mL*d^2/4;
TLx = -1/2*L;
TLyz = 0;

m = mL + M;
Tx = (mL*TLx)/m;
Tyz = 0;
Iyz = ILyz + (TLx - Tx)^2*mL + Tx^2*M;
Ix = 1/2*mL*d^2/4;

T = [Tx;Tyz;Tyz];
I = diag([Ix;Iyz;Iyz]);