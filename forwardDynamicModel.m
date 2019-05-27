function varargout = forwardDynamicModel(varargin)
% vypocita primy dynamicky model manipulatoru 
%   obecne: B(q)ddq + C(q,dq)dq + g(q) = tau - J'(q)h
%     J(q)... kin. jakobian
%     h ... [f;u] ... sila, moment pusobici na konc. eff. (s.s.
%           FN) w.r.t. s.s. F0, h = endEff_force_moment
%   vypocet B(q), tau' = C(q,dq)dq + g(q) + J'(q)h => ddq = B(q)^-1(tau-tau')
%
% varargin = [jointCoords,joint_force_moment,endEff_force_moment,kinPar,dynPar]
%
% jointCoords = [q,dq] ... polohy, rychlosti kloubu
% joint_force_moment = [tau1;tau2;...] ... sily/momenty v kloubech
%   manipulatoru
% endEff_force_moment = [f;u] ... sila, moment pusobici na konc. eff. (s.s.
%   FN) w.r.t. s.s. F0
% kinPar
%   kinPar.DHpar = [d_1,theta_1,a_1,alpha_1; d_2,theta_2,a_2,alpha_2;... ]
%   kinPar.jointType = [joint1, joint2,...] = ['P','R',...]
%   kinPar.Qhome = [q1h; q2h; ...]
% dynPar
%   dynPar.mass = [m1; m2; ...]
%   dynPar.inertiaTensor = [I1, I2, ...] ... w.r.t. s.s. prislusneho 
%       ramene (F1,F2,...) => konstantni (nezavisi na poloze manip.!)
%   dynPar.gravityCenter = [T1, T2, ...] ... w.r.t. s.s. prislusneho 
%       ramene (F1,F2,...) => konstantni (nezavisi na poloze manip.!)
%   dynPar.gravityVector = [0;0;-9.81];
%
% varargout = [tau1;tau2;...]
% [tau1;tau2;...] ... sily/momenty v kloubech manipulatoru

jointCoords = varargin{1};
joint_force_moment = varargin{2};
endEff_force_moment = varargin{3};
kinPar = varargin{4};
dynPar = varargin{5};

DHpar = kinPar.DHpar;
jointType = kinPar.jointType;
Qhome = kinPar.Qhome;

mass = dynPar.mass;
inertiaTensor = dynPar.inertiaTensor;
gravityCenter = dynPar.gravityCenter;
gravityVector = dynPar.gravityVector;

% gravitacni vektor v s.s. F0
g0__0 = gravityVector;

N = size(DHpar,1);

% vypocet matice B(q)
dynPar_zeroGravity = dynPar;
dynPar_zeroGravity.gravityVector = [0;0;0];

ddq_ = zeros(N,1);
for i = 1:N
    ddq_(i) = 1;
    B(:,i) = inverseDynamicModel([jointCoords(:,1),zeros(N,1),ddq_],zeros(6,1),kinPar,dynPar_zeroGravity);
    ddq_(i) = 0;
end

% vypocet tau'
tau_ = inverseDynamicModel([jointCoords(:,1:2),zeros(N,1)],endEff_force_moment,kinPar,dynPar);

% vypocet ddq
ddq = B^-1*(joint_force_moment-tau_);

varargout{1} = ddq;