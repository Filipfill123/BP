function varargout = S2TAV(S,smax,vmax,amax)
% casove optimalni planovani trajektorie s omezenim rychlosti a zrychleni
% varargout = S2TAV(S,smax,vmax,amax)
%
% varargout = [T,A,V], kde T je cas, A je zrychleni, V je rychlost
% odpovidajici bodu na trajektorii v ujete vzdalenosti S
%
% omezeni amax, vmax
% délka trajektorie smax
% 
% REVIZE:
%   2017-01-05: doplnena abs() clenu pod odmocninami (chyba +-1e-16)

% vypocet T, V, A v zavislosti na ujete draze S
if vmax/amax < sqrt(smax/amax)
    S1 = (1/2)*vmax^2/amax;
    S2 = smax-(1/2)*vmax^2/amax;
    S3 = smax;

    if (S <= S1)
        T = sqrt(2)*sqrt(amax*S)/amax;
        V = amax*T;
        A = amax;
    elseif (S>S1)&&(S<S2)
        T = (1/2)*(vmax^2+2*amax*S)/(amax*vmax);
        V = vmax;
        A = 0;
    else
        T = (vmax^2+smax*amax-sqrt(abs(2*vmax^2*smax*amax-2*S*amax*vmax^2)))/(vmax*amax);
        V = -amax*T+amax*(vmax/amax+smax/vmax);
        A = -amax;
    end
else
    S1 = (1/2)*smax;
    S3 = smax;

    if (S <= S1)
        T = sqrt(2)*sqrt(amax*S)/amax;
        V = amax*T;
        A = amax;
    else
        T = -(-2*amax*sqrt(smax/amax)+sqrt(abs(2*smax*amax-2*amax*S)))/amax;
        V = -amax*T+2*amax*sqrt(smax/amax);
        A = -amax;
    end
end
varargout{1} = T;
varargout{2} = A;
varargout{3} = V;