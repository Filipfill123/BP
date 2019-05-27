function varargout = EulerAnglesXYZ2rotMatrix(varargin)

% schema rotace XYZ (EulerAngles = [alpha;beta;gamma])
% prevod [EulerAngles,dEulerAngles,ddEulerangles] -> [rotationMatrix,omega,domega]
% OBSOLETE !!! PREDELAT JAKO EulerAnglesZYX2rotMatrixAndAngularVel.m !!!

EulerAngles = varargin{1};

alpha = EulerAngles(1);
beta = EulerAngles(2);
gamma = EulerAngles(3);

R = [cos(beta)*cos(gamma),-cos(beta)*sin(gamma),sin(beta);sin(alpha)*sin(beta)*cos(gamma)+cos(alpha)*sin(gamma),-sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma),-sin(alpha)*cos(beta);-cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma),cos(alpha)*sin(beta)*sin(gamma)+sin(alpha)*cos(gamma),cos(alpha)*cos(beta)];

varargout{1} = R;

if nargin > 1
    
    dEulerAngles = varargin{2};
    
    H = [1 0 sin(beta);0 cos(alpha) -sin(alpha)*cos(beta);0 sin(alpha) cos(alpha)*cos(beta)];
    omega = H*dEulerAngles;
    
    varargout{2} = omega;
    
    if nargin == 3
        
        dalpha = dEulerAngles(1);
        dbeta = dEulerAngles(2);
        dgamma = dEulerAngles(3);
        
        ddEulerAngles = varargin{3};
        
        dH = [0 0 cos(beta)*dbeta;0 -sin(alpha)*dalpha sin(alpha)*sin(beta)*dbeta-cos(alpha)*cos(beta)*dalpha;0 cos(alpha)*dalpha -cos(alpha)*sin(beta)*dbeta-sin(alpha)*cos(beta)*dalpha];
        
        domega = dH*dEulerAngles + H*ddEulerAngles;
        
        varargout{3} = domega; 
    end
end