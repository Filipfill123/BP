function varargout = kinJacobian(varargin)
% vypocet kinematickeho jakobianu / derivace jakobianu (pro size(jointCoords,2) = 1 / size(jointCoords,2) = 2) 
%
% jacobian = kinJacobian(varargin)
%
% varargin = {jointCoords,DHpar,jointType,baseComp,endEffComp}
%
% varargout = {J,dJ,Jcomp,Jadd} ... pro size(jointCoords,2) = 2 (rychlosti i zrychleni)
% varargout = {J,Jcomp} ... pro size(jointCoords,2) = 1 (pouze rychlosti)
%   J ... kin. jakobian
%   dJ ... derivace J
%   Jcomp, Jadd ... kompenzace: 
%       dX = Jcomp*J*dq
%       ddX = Jcomp*(dJ*dq + J*ddq) + Jadd
%
% jointCoords = [[q1;q2;...],[dq1;dq2;...]]
% DHpar = [d_1,theta1,a_1,alpha_1;d_2,theta_2,a_2,alpha_2; ...]
% jointType = ['R';'P';'R';...]

jointCoords = varargin{1};
DHparam = varargin{2};
jointType = varargin{3};

N = size(DHparam,1);

% nastaveni T0b, TeN - kompenzace polohy zakldny a konc. ef.
if nargin > 3 % pokud kompenzace zadany, jinak bez kompenzaci (= identicke kompenzace)
    baseComp = varargin{4};
    endEffComp = varargin{5};
    
    T0b = [[baseComp(:,2:4),baseComp(:,1)];[0,0,0,1]];
    TeN = [[endEffComp(:,2:4),endEffComp(:,1)];[0,0,0,1]];
else
    T0b = eye(4);
    TeN = eye(4);
end

% vypocet hom. transformacnich matic
q = jointCoords(:,1);

TT{1} = eye(4); % TT{i} = T_{i-1}^0 => TT = {T00,T10,T20,T30,...}
for i = 1:N
    if strcmp(jointType(i),'P')
        DHparam(i,1) = q(i);
        sigma(i) = 1;
    else
        DHparam(i,2) = q(i);
        sigma(i) = 0; 
    end
    sigma_hat(i) = 1 - sigma(i);
    
    TT{i+1} = TT{i}*DH(DHparam(i,:));
end

% vypocet J
for j = 1:N
    z__j_1 = TT{j}(1:3,3);
    r__j_1__N = TT{N+1}(1:3,4) - TT{j}(1:3,4);
    
    if strcmp(jointType(j),'P')
        jP_j = z__j_1;
        jO_j = zeros(3,1);        
    else
        jP_j = cross(z__j_1,r__j_1__N);
        jO_j = z__j_1;
    end
    J(:,j) = [jP_j;jO_j];
end

% vypocet Jcomp (kompenzace)
R0b = T0b(1:3,1:3);
r__N__e__N = TeN(1:3,4);
RN0 = TT{end}(1:3,1:3);

Jcomp = [R0b,-R0b*RN0*scew(r__N__e__N)*RN0';zeros(3,3),R0b];

if size(jointCoords,2) > 1
   % vypocet dJ 

   dq = jointCoords(:,2);
   
   dO_vect = zeros(3,1); % dO_vect = [dO_0,dO_1,...,dO_N]; 
   omega_vect = zeros(3,1); % omega_vect = [omega_0,omega_1,...,omega_N];
   dX_N = J*dq;
   dO_vect(:,N+1) = dX_N(1:3);
   omega_vect(:,N+1) = dX_N(4:6);
   
    for j = 1:N
        z__j_1 = TT{j}(1:3,3);
        r__j_1__j = TT{j+1}(1:3,4) - TT{j}(1:3,4);
        r__j_1__N = TT{N+1}(1:3,4) - TT{j}(1:3,4);
        
        dz__j_1 = cross(omega_vect(:,j),z__j_1);
        dr__j_1__N = dO_vect(:,N+1) - dO_vect(:,j); 
        
        % vypocet sloupcu dJ
        if strcmp(jointType(j),'P')
            djP_j = dz__j_1;
            djO_j = zeros(3,1);        
        else
            djP_j = cross(dz__j_1,r__j_1__N) + cross(z__j_1,dr__j_1__N);
            djO_j = dz__j_1;
        end
        dJ(:,j) = [djP_j;djO_j];
        
        % rekurzivni vypocet rychlosti s.s. (translacni, uhlova)
        omega_vect(:,j+1) = omega_vect(:,j) + z__j_1*sigma_hat(j)*dq(j);
        dO_vect(:,j+1) = dO_vect(:,j) + cross(omega_vect(:,j+1),r__j_1__j) + z__j_1*sigma(j)*dq(j);
    end
    
    % vypocet Jadd(kompenzace)
    omega__N__N = RN0'*omega_vect(:,end);
    
    Jadd = [R0b*RN0*scew(omega__N__N)*scew(omega__N__N)*r__N__e__N;zeros(3,1)];
    
    varargout = {J,dJ,Jcomp,Jadd};
else
    varargout = {J,Jcomp};
end


           



    


