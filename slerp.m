function varargout = slerp(varargin)
% varargin = {Quat,u,mode}
% u = <0,1>
% varargout{1} = [Quat_int,dQuat_du_int,d2Quat_du2_int] ... poloha + derivace dle parametru u
% varargout{2} = u_knot ... uzlovy vektor (pro interpolovane body Quat - rozdeleny dle lin. vzd. mezi int. body)
% mode = 0 ... standardni vypocitany smer rotace 
% mode = 1 ... smer rotace tak, aby byl ujet nejmensi uhel (nejkratsi cesta)

Quat = varargin{1};
u = varargin{2};

if nargin == 3
    mode = varargin{3};
else
    mode = 0;
end

N = size(Quat,2);

% u_knot rozdelen v pomeru ujete vzdalenosti SLERPem
theta_vect(1) = 0;
for i = 1:N-1
    deltaTheta = 2*acos(Quat(:,i)'*Quat(:,i+1));
    if ((mode == 1) && (abs(deltaTheta) > pi))
        Quat(:,i+1) = -Quat(:,i+1);
        deltaTheta = 2*acos(Quat(:,i)'*Quat(:,i+1));
    end
    
    theta_vect = [theta_vect,theta_vect(end) + deltaTheta];
end

u_knot = theta_vect/theta_vect(end);

% ***** Quat, dQuat_dh, d2Quat_dh2 *****
seg = findSpan(u,u_knot);
h = (u - u_knot(seg))/(u_knot(seg+1) - u_knot(seg)); % h = <0,1>

Quat_int = slerp2(Quat(:,seg),Quat(:,seg+1),h);
dQuat_dh_int = quatMult(Quat_int,quatLog(quatMult(quatInv(Quat(:,seg)),Quat(:,seg+1))));
d2Quat_dh2_int = quatMult(Quat_int,quatMult(quatLog(quatMult(quatInv(Quat(:,seg)),Quat(:,seg+1))),quatLog(quatMult(quatInv(Quat(:,seg)),Quat(:,seg+1)))));

dQuat_du_int = dQuat_dh_int/(u_knot(seg+1) - u_knot(seg));
d2Quat_du2_int = d2Quat_dh2_int/(u_knot(seg+1) - u_knot(seg))^2;

varargout{1} = [Quat_int,dQuat_du_int,d2Quat_du2_int];
varargout{2} = u_knot;
varargout{3} = theta_vect;
end


% ********** NESTED FUNCTION **********
function Q = slerp2(Q1,Q2,h)
    % h in <0,1>
    Q = quatMult(Q1,quatExponentiation(quatMult(quatInv(Q1),Q2),h));
end