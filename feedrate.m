function varargout = feedrate(arcLength,u_knot,S_knot,polySegPar)
% arcLength = [S;V;A] ... delka ujeteho oblouku - [poloha;rychlost;zrychleni] [m;m/s;m/s2]
% Pvect,Ovect,Ivect,polyPar ... parametry trajektorie 
% arcLength = [S;V;A] .... profil ujete drahy
% feedType = 1 (aprox) = 2 (presne, num. alg)

S = arcLength(1);
V = arcLength(2);
A = arcLength(3);

% poloha - numericky algoritmus
seg = findSpan(S,S_knot);
numOfSeg = length(S_knot)-1;

if mod(seg,2) == 0
    % poly blend
    [temp,index] = min(abs([S_knot(seg),S_knot(seg+1)] - S));
    u_init = u_knot(seg + index - 1);

    i = seg/2;
    fun_arcLength_integrand = @(u)integrandPoly(u,polySegPar{i});

    err = inf;
    numOfIter = 0;
    u = u_init;
    while err > 1e-4
        if abs(u_knot(seg)-u) < 1e-8
            fun = S_knot(seg);
        else
            fun = S_knot(seg) + quad(fun_arcLength_integrand,u_knot(seg),u);
        end

        err = abs(fun-S);
        dS_du = fun_arcLength_integrand(u);

        delta_u = -dS_du^-1*(fun-S);
        u = u + delta_u;

        numOfIter = numOfIter+1;

        if numOfIter > 20
            arcLength
            error('Prekrocen max. pocet iteraci pri vypoctu feedrate v poly (blending) segmentu')
        end

    end

    % rychlost, zrychleni
    [dS_du,d2S_du2] = fun_arcLength_integrand(u); 
    du_dt = dS_du^-1*V;
    d2u_dt2 = dS_du^-1*(A - d2S_du2*du_dt^2);
else
    % lin. segment
    u = u_knot(seg) + (u_knot(seg+1) - u_knot(seg))*(S - S_knot(seg))/(S_knot(seg+1)-S_knot(seg));

    % rychlost, zrychleni
    du_dt = (u_knot(seg+1) - u_knot(seg))*V/(S_knot(seg+1)-S_knot(seg));
    d2u_dt2 = (u_knot(seg+1) - u_knot(seg))*A/(S_knot(seg+1)-S_knot(seg));

    numOfIter = 0;
end


par_u = [u;du_dt;d2u_dt2]; 

varargout{1} = par_u;
varargout{2} = numOfIter;