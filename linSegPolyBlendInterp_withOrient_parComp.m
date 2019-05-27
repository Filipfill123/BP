function data = linSegPolyBlendInterp_withOrient_parComp(Pvect,Qvect,radiusVect,feedType)
% parametry pro funkci vycisleni interpolovane trajektorie: linSegPolyBlendInterp_withOrient_eval.m v zavislosti na ujete draze
% Pvect = [P1,P2,...], Qvect = [Q1,Q2,...] ... interpolovane body translace, orientace (kvaternion) (feedrate se pocita pouze z translace)
% radiusVect ... vektor radiusu v translaci, kdy zacina blending na linearnich segmentech (vzd. okamziku blendovani od vrcholu Pi)
%   Interpolace orientace v lin. segmentech = SLERP (od Ovect do Ivect - pocatek konec blendovani v translaci) + prima interpolace prvku kvaternionu (ne zcela dobre) v blending segmentech kubickym SPLINEM (od Ivect do Ovect)
% res ... rozliseni ujete drahy v [m] 
% feedType  = 1 (aproximativni vypocet feedrate = fluktuace tecne rychlosti, zpomaleni v blending segmentech - jednoduchy ALG., aprox. prirozena parametrizace)
%           = 2 (dodrzeni predepsaneho profilu rychlosti - presna zavislost na ujete draze, numer. ALG, presna prirozena parametrizace)
% 

numOfPoints = size(Pvect,2);

% POZOR: osetreni stejnych koincidencnich bodu (spravne by to bylo osetrit dale v algoritmu! ZDE: pokud jsou dva po sobe jdouci koincidencni body shodne, jeden se perturbuje malou hodnotou)
for i = 1:size(Pvect,2)-1
    % translace
    if norm(Pvect(:,i)-Pvect(:,i+1)) < 1e-8
        disp('******** IDENTICKE ZADANE TRANSLACE !!! **********')
        disp(['Pvect(:,',num2str(i),') = Pvect(:,',num2str(i+1),') = [',num2str(Pvect(:,i)'),']'])
        % mala perturbace
        Pvect(:,i+1) = Pvect(:,i+1) + 1e-5*rand(3,1);
        disp([' => (perturbace) Pvect(:,',num2str(i+1),') = [',num2str(Pvect(:,i+1)'),']'])
        disp('******************')
    end
    % orientace
    dQuat = quatMult(Qvect(:,i),quatInv(Qvect(:,i+1)));
    if abs(abs(dQuat(1))-1) < 1e-8
        disp('******** IDENTICKE ZADANE ORIENTACE !!! **********')
        disp(['Qvect(:,',num2str(i),') = Qvect(:,',num2str(i+1),') = [',num2str(Qvect(:,i)'),']'])
        % mala perturbace
        Qvect(:,i+1) = Qvect(:,i+1) + 1e-5*rand(4,1);
        Qvect(:,i+1) = Qvect(:,i+1)/norm(Qvect(:,i+1));
        disp([' => (perturbace) Qvect(:,',num2str(i+1),') = [',num2str(Qvect(:,i+1)'),']'])
        disp('******************')
    end
    
end

% osetreni spatne nastavenych radiusu
for i = 2:length(radiusVect)-1
    if norm(Pvect(:,i+1)-Pvect(:,i)) <= radiusVect(i-1) + radiusVect(i)
        error(['Radius ',num2str(i-1),' a/nebo ',num2str(i),' je prilis velky.'])
    end
end
if norm(Pvect(:,2)-Pvect(:,1)) <= radiusVect(1);
    error(['Prvni radius je prilis velky.'])
end
if norm(Pvect(:,end)-Pvect(:,end-1)) <= radiusVect(end);
    error(['Posledni radius je prilis velky.'])
end

% pocatky a konce blendovani
% translace
for i = 2:numOfPoints-1
    % translace
    Ovect(:,i-1) = (Pvect(:,i+1)-Pvect(:,i))/norm(Pvect(:,i+1)-Pvect(:,i))*radiusVect(i-1) + Pvect(:,i);
    Ivect(:,i-1) = -(Pvect(:,i)-Pvect(:,i-1))/norm(Pvect(:,i)-Pvect(:,i-1))*radiusVect(i-1) + Pvect(:,i);
end

% orientace (parametr theta pro zacatek a konec blendovani (pro Ivect, Ovect) vyplyva z SLERP inter. mezi zadanymi kvaterniony v koincidencnich bodech)
% interpolace v lin. segmentech = SLERP (od Ovect do Ivect), prima interpolace prvku kvaternionu (ne zcela dobre) v blending segmentech kubickym SPLINEM (od Ivect do Ovect)
for i = 1:numOfPoints-1
    deltaTheta(:,i) = 2*acos(Qvect(:,i)'*Qvect(:,i+1));
    if (abs(deltaTheta(:,i)) > pi) % prejezd "kratsi" stranou
        Qvect(:,i+1) = -Qvect(:,i+1);
        deltaTheta(:,i) = 2*acos(Qvect(:,i)'*Qvect(:,i+1));
    end
end

for i = 1:numOfPoints-2
    deltaThetaIvect(:,i) = (1-norm(Pvect(:,i+1)-Ivect(:,i))/norm(Pvect(:,i+1)-Pvect(:,i)))*deltaTheta(:,i); % delka o kterou ujede theta na linearnim segmentu
    deltaThetaOvect(:,i) = (norm(Pvect(:,i+1)-Ovect(:,i))/norm(Pvect(:,i+2)-Pvect(:,i+1)))*deltaTheta(:,i+1);
end

% delky segmentu - APROXIMACE
arcLength_polySeg_approx = polySegLength_approx(Pvect,Ovect,Ivect);
arcLength_linSeg = linSegLength(Pvect,Ovect,Ivect);

% APROXIMACE: ujeta draha + uzlovy vektor u = <0,1>
S_knot = [0];
for i = 1:numOfPoints-2
    S_knot = [S_knot, S_knot(end) + arcLength_linSeg(i)];
    S_knot = [S_knot,S_knot(end) + arcLength_polySeg_approx(i)];
end
S_knot_approx = [S_knot,S_knot(end) + arcLength_linSeg(end)];
        
S_max_approx = S_knot_approx(end);
u_knot = S_knot_approx/S_max_approx;

[linSegPar_trans,linSegPar_orient] = linSegParam(u_knot,Pvect,Ovect,Ivect,deltaThetaOvect,deltaThetaIvect,deltaTheta,Qvect);
[polySegPar_trans,polySegPar_orient] = polySegParam(u_knot,linSegPar_trans,linSegPar_orient,Qvect,deltaTheta);

% vypocet feedrate (dle TRANSLACE!!!!!) + generovani trajektorie
switch feedType
    case 1 % aproximace
        S_knot = S_knot_approx;

    case 2 % presne
        
        % skutecna delka polynomialnich segmentu
        S_knot = [0];
        arcLength_polySeg = polySegLength(u_knot,polySegPar_trans);
        for i = 1:numOfPoints-2
            S_knot = [S_knot, S_knot(end) + arcLength_linSeg(i)];
            S_knot = [S_knot,S_knot(end) + arcLength_polySeg(i)];
        end
        S_knot = [S_knot,S_knot(end) + arcLength_linSeg(end)];
        S_max = S_knot(end);

    otherwise
        error('Spatny feedType')
end

data.feedType = feedType; % typ vypoctu feedrate
data.u_knot = u_knot; % uzlovy vektor
data.S_knot = S_knot; % ujeta draha (v translaci) v uzlovem vektoru (lin a poly segmenty)
data.deltaTheta = deltaTheta; % ujeta draha v orientaci mezi koincidencnimi body
data.linSegPar_trans = linSegPar_trans; % parametry generovane trajektorie
data.linSegPar_orient = linSegPar_orient; % parametry generovane trajektorie
data.polySegPar_trans = polySegPar_trans; % parametry generovane trajektorie
data.polySegPar_orient = polySegPar_orient; % parametry generovane trajektorie
data.Ivect = Ivect;
data.Ovect = Ovect;