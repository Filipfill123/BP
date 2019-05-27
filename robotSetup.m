function Links = robotSetup(varargin)

kinPar = varargin{1};
dynPar = varargin{2};

DHpar = kinPar.DHpar;
jointType = kinPar.jointType;
Qhome = kinPar.Qhome;
endEffComp = kinPar.endEffComp;
baseComp = kinPar.baseComp;

mass = dynPar.mass;
inertiaTensor = dynPar.inertiaTensor;
gravityCenter = dynPar.gravityCenter;
gravityVector = dynPar.gravityVector;

N = size(DHpar,1);

% kompenzace polohy konc. ef.
Oen = endEffComp(:,1);
Ren = endEffComp(:,2:4);

% kompenzace polohy zakladny
O0b = baseComp(:,1);
R0b = baseComp(:,2:4);

% kompenzace polohy zakladny
TT{1} = [[R0b,O0b];[0,0,0,1]]; % T0b
Links{1}.A = [0;0;0];
Links{1}.B = TT{1}(1:3,4);
Links{1}.R = TT{1}(1:3,1:3);
Links{1}.R0b = R0b;
Links{1}.O0b = O0b;
Links{1}.gravityVector = gravityVector; 
Links{1}.numOfLinks = N;
    
TT_b{1} = TT{1}; % TT_b = {T0b, T1b, T2b, ...}
for i = 1:N; 
    TT{i+1} = DH(DHpar(i,:)); % TT = {T0b,T10,T21, ...}
    TT_b{i+1} = TT_b{i}*TT{i+1};

    % Links Links = {baseComp=Link0,Link1,Link2,...LinkN, LinkN+1=endEffComp}
    Links{i+1}.A = TT_b{i}(1:3,4);
    Links{i+1}.B = TT_b{i+1}(1:3,4);
    Links{i+1}.R = TT_b{i+1}(1:3,1:3);
    Links{i+1}.qhome = Qhome(i);
    Links{i+1}.mass = mass(i);
    Links{i+1}.inertiaTensor = inertiaTensor(:,(3*(i-1)+1):(3*(i-1)+3));
    GC = TT_b{i+1}*[gravityCenter(:,i);1];
    Links{i+1}.gravityCenter = GC(1:3);
    Links{i+1}.jointAxis = TT_b{i}(1:3,3);

    if strcmp(jointType(i),'R') % sigma = 0 for R joint, sigma = 1 for P joint
        Links{i+1}.jointType_sigma = 0;
    else
        Links{i+1}.jointType_sigma = 1;
    end   
end

% kompenzace polohy konc. efektoru
TT{N+2} = [[Ren,Oen];[0,0,0,1]];
TT_b{N+2} = TT_b{N+1}*TT{N+2};
Links{N+2}.A = TT_b{N+1}(1:3,4);
Links{N+2}.B = TT_b{N+2}(1:3,4);
Links{N+2}.R = TT_b{i+2}(1:3,1:3);
Links{N+2}.Ren = Ren;
Links{N+2}.Oen = Oen;

