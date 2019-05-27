function varargout = linSegPolyBlend_withOrient_feedrate(S,u_knot,S_knot,polySegPar)
% S ... delka ujeteho oblouku - [m]
% u_knot,S_knot,polySegPar ... parametry trajektorie 
% feedType = 1 (aprox) = 2 (presne, num. alg)

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

    % rychlost, zrychleni dle promenne S
    [dS_du,d2S_du2] = fun_arcLength_integrand(u); 
    du_dS = dS_du^-1;
    d2u_dS2 = -dS_du^-3*d2S_du2;
else
    % lin. segment
    u = u_knot(seg) + (u_knot(seg+1) - u_knot(seg))*(S - S_knot(seg))/(S_knot(seg+1)-S_knot(seg));

    % rychlost, zrychleni dle promenne S
    du_dS = (u_knot(seg+1) - u_knot(seg))/(S_knot(seg+1)-S_knot(seg));
    d2u_dS2 = 0;

    numOfIter = 0;
end


par_u = [u;du_dS;d2u_dS2]; 

varargout{1} = par_u;
varargout{2} = numOfIter;